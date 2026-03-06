#ifndef SOA_BUFFER_HPP
#define SOA_BUFFER_HPP
// =============================================================================
//  SoABuffer.hpp
//  Automatic Structure-of-Arrays buffer from any C++ POD struct,
//  with ROOT TTree branch registration and clear support.
//
//  Requirements:
//    - C++17 or later  (GCC >= 12 supported)
//    - Boost >= 1.75   (Boost.PFR, header-only)
//    - ROOT >= 6.x     (for TTree / TBranch)
//
//  Field-name reflection:
//    In C++20 mode (Boost >= 1.80) real struct field names are used for
//    branch names (e.g. "trk_px").  In C++17 mode the fallback scheme
//    "field_0", "field_1", ... is used automatically.
//
//  Compile (C++17):
//    g++ -std=c++17 main.cpp $(root-config --cflags --libs) -I/path/to/boost -o soa_demo
//  Compile (C++20, real field names):
//    g++ -std=c++20 main.cpp $(root-config --cflags --libs) -I/path/to/boost -o soa_demo
// =============================================================================

#include <boost/pfr.hpp>
#include <TTree.h>

#include <array>
#include <iostream>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>
#include <stdexcept>

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------
namespace soa_detail {

// Map a tuple<T0, T1, ...> → tuple<vector<T0>, vector<T1>, ...>
template<typename Tuple>
struct TupleToVectors;

template<typename... Ts>
struct TupleToVectors<std::tuple<Ts...>> {
    using type = std::tuple<std::vector<Ts>...>;
};

// Map a tuple<T0, T1, ...> → tuple<vector<T0>*, vector<T1>*, ...>
// Used to provide the T** that ROOT's SetBranchAddress requires.
template<typename Tuple>
struct TupleToVectorPtrs;

template<typename... Ts>
struct TupleToVectorPtrs<std::tuple<Ts...>> {
    using type = std::tuple<std::vector<Ts>*...>;
};

// Build a runtime array of field name strings.
// C++20 + Boost >= 1.80 : real field names via pfr::names_as_array.
// C++17 fallback        : "field_0", "field_1", ...
template<typename Struct, std::size_t... Is>
std::array<std::string, sizeof...(Is)>
make_field_names_impl(std::index_sequence<Is...>) {
#if __cplusplus >= 202002L
    constexpr auto pfr_names = boost::pfr::names_as_array<Struct>();
    return { std::string(pfr_names[Is])... };
#else
    return { ("field_" + std::to_string(Is))... };
#endif
}

template<typename Struct>
std::array<std::string, boost::pfr::tuple_size_v<Struct>>
get_field_names() {
    return make_field_names_impl<Struct>(
        std::make_index_sequence<boost::pfr::tuple_size_v<Struct>>{});
}

} // namespace soa_detail

// ---------------------------------------------------------------------------
//  SoABuffer<Struct>
//  ---
//  Wraps one std::vector<T> per field of Struct, reflecting the layout
//  automatically through Boost.PFR.  Provides:
//    push_back(const Struct&)   – append one element
//    get(size_t i)              – reconstruct an AoS element
//    size()                     – number of stored elements
//    clear()                    – empty all vectors (keeps capacity)
//    make_branches(TTree*, prefix)  – register all vectors as TTree branches
//    set_branch_addresses(TTree*, prefix) – re-point addresses on an existing tree
// ---------------------------------------------------------------------------
template<typename Struct>
class SoABuffer {
    static_assert(std::is_trivially_copyable_v<Struct>,
                  "SoABuffer requires a POD / trivially-copyable struct");

public:
    // Number of fields in Struct
    static constexpr std::size_t kNFields = boost::pfr::tuple_size_v<Struct>;

    // Tuple-of-vectors type that mirrors the struct layout
    using FieldTuple  = decltype(boost::pfr::structure_to_tuple(std::declval<Struct>()));
    using ArraysTuple = typename soa_detail::TupleToVectors<FieldTuple>::type;
    // Tuple of raw pointers to each vector — passed to ROOT's SetBranchAddress (needs T**)
    using PtrsTuple   = typename soa_detail::TupleToVectorPtrs<FieldTuple>::type;

    // ------------------------------------------------------------------
    // Construction
    // ------------------------------------------------------------------
    SoABuffer() : ptrs_(make_ptrs(std::make_index_sequence<kNFields>{})) {}

    // Copying would leave ptrs_ pointing at the source's arrays_ — disallow.
    SoABuffer(const SoABuffer&)            = delete;
    SoABuffer& operator=(const SoABuffer&) = delete;

    // Move is safe: re-initialise ptrs_ to point at the new arrays_.
    SoABuffer(SoABuffer&& o) noexcept
        : arrays_(std::move(o.arrays_))
        , ptrs_(make_ptrs(std::make_index_sequence<kNFields>{})) {}

    SoABuffer& operator=(SoABuffer&& o) noexcept {
        if (this != &o) {
            arrays_ = std::move(o.arrays_);
            ptrs_   = make_ptrs(std::make_index_sequence<kNFields>{});
        }
        return *this;
    }

    /// Reserve memory for all vectors at once
    void reserve(std::size_t n) {
        reserve_impl(n, std::make_index_sequence<kNFields>{});
    }

    // ------------------------------------------------------------------
    // Element access
    // ------------------------------------------------------------------

    /// Append a struct as a new row
    void push_back(const Struct& s) {
        push_impl(s, std::make_index_sequence<kNFields>{});
    }

    /// Reconstruct a struct from row i
    Struct get(std::size_t i) const {
        return get_impl(i, std::make_index_sequence<kNFields>{});
    }

    /// Number of rows stored
    std::size_t size() const {
        return std::get<0>(arrays_).size();
    }

    /// Direct access to the i-th column vector (type-safe via index)
    template<std::size_t I>
    auto& column() { return std::get<I>(arrays_); }

    template<std::size_t I>
    const auto& column() const { return std::get<I>(arrays_); }

    // ------------------------------------------------------------------
    // Clear
    // ------------------------------------------------------------------

    /// Empty all column vectors (capacity is preserved)
    void clear() {
        clear_impl(std::make_index_sequence<kNFields>{});
    }

    // ------------------------------------------------------------------
    // ROOT TTree interface
    // ------------------------------------------------------------------

    /// Create one STL-vector branch per field on *tree*.
    /// The branch name is  prefix + field_name  (e.g. "trk_x", "trk_y", …).
    /// Call once after TTree construction, before the event loop.
    void make_branches(TTree* tree, const std::string& prefix = "") {
        if (!tree)
            throw std::invalid_argument("SoABuffer::make_branches: null TTree pointer");
        make_branches_impl(tree, prefix, std::make_index_sequence<kNFields>{});
    }

    /// Re-point branch addresses to this buffer's vectors.
    /// Use when reading back from an existing file, or after the buffer
    /// has been moved in memory.
    void set_branch_addresses(TTree* tree, const std::string& prefix = "") {
        if (!tree)
            throw std::invalid_argument("SoABuffer::set_branch_addresses: null TTree pointer");
        set_addresses_impl(tree, prefix, std::make_index_sequence<kNFields>{});
    }

    // ------------------------------------------------------------------
    // Utilities
    // ------------------------------------------------------------------

    /// Field names — real names in C++20 mode, "field_N" in C++17 mode.
    static std::array<std::string, kNFields> field_names() {
        return soa_detail::get_field_names<Struct>();
    }

    /// Print a short summary of stored columns to stdout
    void print_summary(std::ostream& os = std::cout) const {
        auto names = soa_detail::get_field_names<Struct>();
        os << "SoABuffer<" << typeid(Struct).name()
           << ">  rows=" << size()
           << "  fields=" << kNFields << '\n';
        print_impl(os, names, std::make_index_sequence<kNFields>{});
    }

private:
    ArraysTuple arrays_;
    PtrsTuple   ptrs_;   // each element points at the corresponding vector in arrays_

    template<std::size_t... Is>
    PtrsTuple make_ptrs(std::index_sequence<Is...>) {
        return PtrsTuple{ &std::get<Is>(arrays_)... };
    }

    // ---- push_back -------------------------------------------------------
    template<std::size_t... Is>
    void push_impl(const Struct& s, std::index_sequence<Is...>) {
        (std::get<Is>(arrays_).push_back(boost::pfr::get<Is>(s)), ...);
    }

    // ---- get -------------------------------------------------------------
    template<std::size_t... Is>
    Struct get_impl(std::size_t i, std::index_sequence<Is...>) const {
        Struct s{};
        ((boost::pfr::get<Is>(s) = std::get<Is>(arrays_)[i]), ...);
        return s;
    }

    // ---- clear -----------------------------------------------------------
    template<std::size_t... Is>
    void clear_impl(std::index_sequence<Is...>) {
        (std::get<Is>(arrays_).clear(), ...);
    }

    // ---- reserve ---------------------------------------------------------
    template<std::size_t... Is>
    void reserve_impl(std::size_t n, std::index_sequence<Is...>) {
        (std::get<Is>(arrays_).reserve(n), ...);
    }

    // ---- make_branches ---------------------------------------------------
    // ROOT TTree::Branch for std::vector<T> takes a pointer-to-pointer:
    //   tree->Branch("name", &vec_ptr)
    // where vec_ptr is a std::vector<T>*.
    template<std::size_t... Is>
    void make_branches_impl(TTree* tree,
                            const std::string& prefix,
                            std::index_sequence<Is...>) {
        auto names = soa_detail::get_field_names<Struct>();
        (tree->Branch(
            (prefix + names[Is]).c_str(),
            &std::get<Is>(arrays_)          // ROOT takes std::vector<T>* directly
        ), ...);
    }

    // ---- set_branch_addresses --------------------------------------------
    // SetBranchAddress for STL-vector branches requires T** (pointer-to-pointer).
    // ptrs_[I] is a std::vector<T>* pointing at arrays_[I]; we pass &ptrs_[I].
    template<std::size_t... Is>
    void set_addresses_impl(TTree* tree,
                            const std::string& prefix,
                            std::index_sequence<Is...>) {
        auto names = soa_detail::get_field_names<Struct>();
        (tree->SetBranchAddress(
            (prefix + names[Is]).c_str(),
            &std::get<Is>(ptrs_)            // T** — what ROOT requires for vector branches
        ), ...);
    }

    // ---- print -----------------------------------------------------------
    template<std::size_t... Is>
    void print_impl(std::ostream& os,
                    const std::array<std::string, kNFields>& names,
                    std::index_sequence<Is...>) const {
        ((os << "  [" << Is << "] " << names[Is]
             << "  size=" << std::get<Is>(arrays_).size() << '\n'), ...);
    }
};

#endif // SOA_BUFFER_HPP