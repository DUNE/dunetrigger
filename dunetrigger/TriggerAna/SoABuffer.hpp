#ifndef SOA_BUFFER_HPP
#define SOA_BUFFER_HPP
// =============================================================================
//  SoABuffer.hpp
//  Automatic Structure-of-Arrays buffer from a C++ struct,
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

} // namespace soa_detail

// ---------------------------------------------------------------------------
//  SoaFieldNames<Struct> — specialisable trait for field name registration.
//
//  C++20: names are always available automatically via Boost.PFR; no action
//         needed from the user.
//
//  C++17: field name reflection does not exist in the language.  Users must
//         explicitly specialise this trait for each struct, or use the
//         convenience macro REGISTER_SOA_FIELD_NAMES.  If no specialisation
//         is provided the fallback "field_0", "field_1", ... names are used
//         and a compile-time warning is emitted via static_assert.
//
//  Usage (C++17, outside any namespace):
//    REGISTER_SOA_FIELD_NAMES(Track, x, y, z, px, py, pz, chi2, n_hits, pdg_id, is_primary)
// ---------------------------------------------------------------------------
template<typename Struct>
struct SoaFieldNames {
    static constexpr bool registered = false;

    static std::array<std::string, boost::pfr::tuple_size_v<Struct>> get() {
        // Fallback: "field_0", "field_1", ...
        return get_impl(std::make_index_sequence<boost::pfr::tuple_size_v<Struct>>{});
    }
private:
    template<std::size_t... Is>
    static std::array<std::string, sizeof...(Is)> get_impl(std::index_sequence<Is...>) {
        return { ("field_" + std::to_string(Is))... };
    }
};

// ---------------------------------------------------------------------------
//  REGISTER_SOA_FIELD_NAMES(StructType, field1, field2, ...)
//  Specialises SoaFieldNames for StructType with the given field names.
//  Place at namespace scope (outside any class or function).
//
//  Uses Boost.Preprocessor to stringify each field name individually:
//
//    BOOST_PP_VARIADIC_TO_SEQ(x, y, z)  →  (x)(y)(z)      [a PP sequence]
//    BOOST_PP_SEQ_ENUM(                  →  "x", "y", "z"  [comma-separated]
//        BOOST_PP_SEQ_TRANSFORM(
//            SOA_PP_STRINGIFY_OP, _, seq))
//
//  BOOST_PP_SEQ_TRANSFORM applies SOA_PP_STRINGIFY_OP(_, _, elem) to every
//  element of the sequence, which expands to BOOST_PP_STRINGIZE(elem) i.e.
//  the quoted token.  BOOST_PP_SEQ_ENUM then joins them with commas so they
//  form a valid brace-initialiser for the const char* array.
//
//  Field limit: BOOST_PP_LIMIT_SEQ (default 256).
// ---------------------------------------------------------------------------
#include <boost/preprocessor/variadic/to_seq.hpp>
#include <boost/preprocessor/seq/transform.hpp>
#include <boost/preprocessor/seq/enum.hpp>
#include <boost/preprocessor/stringize.hpp>

// Callback for BOOST_PP_SEQ_TRANSFORM: (data, elem) → "elem"
// The macro receives (r, data, elem); data is unused (_).
#define SOA_PP_STRINGIFY_OP(r, _, elem) BOOST_PP_STRINGIZE(elem)

// Produce a comma-separated list of quoted strings from a variadic list:
//   SOA_PP_STRINGIFY_EACH(x, y, z)  →  "x", "y", "z"
#define SOA_PP_STRINGIFY_EACH(...)                          \
    BOOST_PP_SEQ_ENUM(                                      \
        BOOST_PP_SEQ_TRANSFORM(                             \
            SOA_PP_STRINGIFY_OP, _,                         \
            BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__)))

#define REGISTER_SOA_FIELD_NAMES(StructType, ...)                               \
template<>                                                                       \
struct SoaFieldNames<StructType> {                                               \
    static constexpr bool registered = true;                                     \
    static constexpr std::size_t kN = boost::pfr::tuple_size_v<StructType>;     \
    static std::array<std::string, kN> get() {                                  \
        static const char* const names[] = {                                    \
            SOA_PP_STRINGIFY_EACH(__VA_ARGS__)                                   \
        };                                                                       \
        constexpr std::size_t kProvided =                                        \
            sizeof(names) / sizeof(names[0]);                                    \
        static_assert(kProvided == kN,                                           \
            "REGISTER_SOA_FIELD_NAMES: name count does not match "              \
            "the number of fields in " #StructType ". "                         \
            "Check that every field has exactly one name.");                     \
        return get_impl(std::make_index_sequence<kN>{}, names);                  \
    }                                                                            \
private:                                                                         \
    template<std::size_t... Is>                                                  \
    static std::array<std::string, sizeof...(Is)>                               \
    get_impl(std::index_sequence<Is...>, const char* const* n) {                \
        return { std::string(n[Is])... };                                        \
    }                                                                            \
};

namespace soa_detail {

// Dispatch: C++20 uses PFR directly; C++17 goes through the trait.
#if __cplusplus >= 202002L
// Forward declaration — defined below get_field_names to avoid lookup issues
// across compilers.
template<typename Struct, typename NamesArray, std::size_t... Is>
std::array<std::string, sizeof...(Is)>
get_names_impl(const NamesArray& pfr_names, std::index_sequence<Is...>);
#endif

template<typename Struct>
std::array<std::string, boost::pfr::tuple_size_v<Struct>>
get_field_names() {
#if __cplusplus >= 202002L
    constexpr auto pfr_names = boost::pfr::names_as_array<Struct>();
    return get_names_impl<Struct>(pfr_names,
        std::make_index_sequence<boost::pfr::tuple_size_v<Struct>>{});
#else
    static_assert(SoaFieldNames<Struct>::registered,
        "SoABuffer (C++17): field names not registered for this struct. "
        "Use REGISTER_SOA_FIELD_NAMES(StructType, field1, field2, ...) "
        "at namespace scope, or compile with -std=c++20.");
    return SoaFieldNames<Struct>::get();
#endif
}

#if __cplusplus >= 202002L
template<typename Struct, typename NamesArray, std::size_t... Is>
std::array<std::string, sizeof...(Is)>
get_names_impl(const NamesArray& pfr_names, std::index_sequence<Is...>) {
    return { std::string(pfr_names[Is])... };
}
#endif

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
    static_assert(std::is_default_constructible_v<Struct>,
                  "SoABuffer requires a default-constructible struct");
    static_assert(std::is_copy_constructible_v<Struct>,
                  "SoABuffer requires a copy-constructible struct");
    static_assert(boost::pfr::tuple_size_v<Struct> > 0,
                  "SoABuffer does not support empty structs");

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

    /// Reconstruct a struct from row i. Throws std::out_of_range if i >= size().
    Struct get(std::size_t i) const {
        if (i >= size())
            throw std::out_of_range(
                "SoABuffer::get: index " + std::to_string(i) +
                " out of range (size=" + std::to_string(size()) + ")");
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

    /// Create one STL-vector branch per field on tree.
    /// The branch name is  prefix + field_name  (e.g. "trk_x", "trk_y", …).
    /// Call once after TTree construction, before the event loop.
    void make_branches(TTree& tree, const std::string& prefix = "") {
        make_branches_impl(&tree, prefix, std::make_index_sequence<kNFields>{});
    }

    /// Re-point branch addresses to this buffer's vectors.
    /// Use when reading back from an existing file, or after the buffer
    /// has been moved in memory.
    void set_branch_addresses(TTree& tree, const std::string& prefix = "") {
        set_addresses_impl(&tree, prefix, std::make_index_sequence<kNFields>{});
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
    // Return value is checked: ROOT returns <0 on failure (name/type mismatch).
    template<std::size_t... Is>
    void set_addresses_impl(TTree* tree,
                            const std::string& prefix,
                            std::index_sequence<Is...>) {
        auto names = soa_detail::get_field_names<Struct>();
        (check_set_address(
            tree->SetBranchAddress(
                (prefix + names[Is]).c_str(),
                &std::get<Is>(ptrs_)        // T** — what ROOT requires for vector branches
            ),
            prefix + names[Is]
        ), ...);
    }

    static void check_set_address(Int_t status, const std::string& branch_name) {
        // ROOT returns kMissingBranch (-5) or other negative codes on failure
        if (status < 0)
            throw std::runtime_error(
                "SoABuffer::set_branch_addresses: failed to bind branch \"" +
                branch_name + "\" (ROOT error code " + std::to_string(status) + ")");
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

// ---------------------------------------------------------------------------
//  SoAWriter<Struct>
//  ---
//  Convenience wrapper around SoABuffer<Struct> with an internal staging row.
//  Typical usage:
//
//    SoAWriter<Track> writer;
//    writer.make_branches(tree, "trk_");
//
//    for (auto& raw : source) {
//        writer->x         = raw.x;      // fill staging row
//        writer->px        = raw.px;
//        writer.push_back();             // commit row → SoA buffer
//    }
//    tree.Fill();
//    writer.clear();                     // reset buffer and staging row
//
//  Methods:
//    push_back()          – append current value of `row` into the SoA buffer
//    commit_and_reset()   – append current `row`, then reset it
//    clear()              – clear the SoA buffer AND zero-initialise `row`
//    reset_row()          – zero-initialise `row` only (buffer unchanged)
//    buffer()             – access the underlying SoABuffer (e.g. for size(),
//                           get(i), set_branch_addresses for read-back)
//    make_branches(...)   – forwarded to SoABuffer
// ---------------------------------------------------------------------------
template<typename Struct>
class SoAWriter {
public:
    /// Access the staging row fields via `writer->field`.
    Struct* operator->() noexcept { return &row_; }
    const Struct* operator->() const noexcept { return &row_; }
    Struct& row() noexcept { return row_; }
    const Struct& row() const noexcept { return row_; }

    // ------------------------------------------------------------------
    // Construction
    // ------------------------------------------------------------------
    SoAWriter() = default;

    explicit SoAWriter(std::size_t reserve_n) {
        buffer_.reserve(reserve_n);
    }

    // ------------------------------------------------------------------
    // Core operations
    // ------------------------------------------------------------------

    /// Copy the current value of `row` into the SoA buffer.
    void push_back() {
        buffer_.push_back(row_);
    }

    /// Copy current row into the buffer and reset staging row.
    void commit_and_reset() {
        push_back();
        reset_row();
    }

    /// Empty the SoA buffer and zero-initialise the staging row.
    void clear() {
        buffer_.clear();
        reset_row();
    }

    /// Zero-initialise the staging row without touching the buffer.
    void reset_row() {
        row_ = Struct{};
    }

    // ------------------------------------------------------------------
    // Buffer access
    // ------------------------------------------------------------------

    [[nodiscard]] SoABuffer<Struct>& buffer() noexcept { return buffer_; }
    [[nodiscard]] const SoABuffer<Struct>& buffer() const noexcept { return buffer_; }

    /// Current SoA buffer row count (number of committed rows).
    [[nodiscard]] std::size_t row_count() const noexcept { return buffer_.size(); }

    /// Shorthand for the number of committed rows.
    [[nodiscard]] std::size_t size() const noexcept { return buffer_.size(); }

    // ------------------------------------------------------------------
    // ROOT TTree forwarding
    // ------------------------------------------------------------------

    void make_branches(TTree& tree, const std::string& prefix = "") {
        buffer_.make_branches(tree, prefix);
    }

    void print_summary(std::ostream& os = std::cout) const {
        buffer_.print_summary(os);
    }

private:
    Struct row_{};
    SoABuffer<Struct> buffer_;
};

#endif // SOA_BUFFER_HPP
