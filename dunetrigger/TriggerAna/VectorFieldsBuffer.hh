#ifndef VECTOR_FIELDS_BUFFER_HH
#define VECTOR_FIELDS_BUFFER_HH
// =============================================================================
//  VectorFieldsBuffer.hh
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
//    branch names (e.g. "trk_px").  In C++17 mode use
//    REGISTER_FIELD_NAMES(StructType, field1, ...) at namespace scope.
//
//  Compile (C++17):
//    g++ -std=c++17 main.cpp $(root-config --cflags --libs) -I/path/to/boost -o soa_demo
//  Compile (C++20, real field names):
//    g++ -std=c++20 main.cpp $(root-config --cflags --libs) -I/path/to/boost -o soa_demo
// =============================================================================

#include "FieldNames.hh"

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

// Map a tuple<T0, T1, ...> -> tuple<vector<T0>, vector<T1>, ...>
template<typename Tuple>
struct TupleToVectors;

template<typename... Ts>
struct TupleToVectors<std::tuple<Ts...>> {
    using type = std::tuple<std::vector<Ts>...>;
};

// Map a tuple<T0, T1, ...> -> tuple<vector<T0>*, vector<T1>*, ...>
// Used to provide the T** that ROOT's SetBranchAddress requires.
template<typename Tuple>
struct TupleToVectorPtrs;

template<typename... Ts>
struct TupleToVectorPtrs<std::tuple<Ts...>> {
    using type = std::tuple<std::vector<Ts>*...>;
};

} // namespace soa_detail

// ---------------------------------------------------------------------------
//  VectorFieldsBuffer<Struct>
//  ---
//  Wraps one std::vector<T> per field of Struct, reflecting the layout
//  automatically through Boost.PFR.
//
//  Write path -- fill the public staging row, then commit it:
//    buf->field = value;    // or buf.row.field = value
//    buf.push_back();       // commit row to SoA storage
//    tree.Fill();
//    buf.clear();           // empty storage and reset row
//
//  Read path -- attach to an existing tree, then iterate:
//    buf.set_branch_addresses(tree, "prefix_");
//    for (Long64_t e = 0; e < tree.GetEntries(); ++e) {
//        tree.GetEntry(e);
//        buf.get(i);        // reconstruct one AoS element
//    }
//
//  Methods:
//    row                             -- public staging struct (read/write)
//    push_back()                     -- commit row to storage (gated)
//    push_back(const Struct&)        -- append directly (always)
//    commit_and_reset()              -- push_back() + reset_row() (gated)
//    reset_row()                     -- zero-initialise row (always)
//    get(size_t i)                   -- reconstruct an AoS element
//    size()                          -- number of stored rows
//    clear()                         -- empty storage and reset row (gated)
//    reserve(n)                      -- pre-allocate column vectors
//    make_branches(TTree&, prefix)   -- register TTree branches (gated)
//    set_branch_addresses(TTree&, prefix) -- attach to existing branches
//    enable(bool) / operator bool()  -- enable/disable write operations
// ---------------------------------------------------------------------------
template<typename Struct>
class VectorFieldsBuffer {
    static_assert(std::is_default_constructible_v<Struct>,
                  "VectorFieldsBuffer requires a default-constructible struct");
    static_assert(std::is_copy_constructible_v<Struct>,
                  "VectorFieldsBuffer requires a copy-constructible struct");
    static_assert(boost::pfr::tuple_size_v<Struct> > 0,
                  "VectorFieldsBuffer does not support empty structs");

public:
    static constexpr std::size_t kNFields = boost::pfr::tuple_size_v<Struct>;

    using FieldTuple  = decltype(boost::pfr::structure_to_tuple(std::declval<Struct>()));
    using ArraysTuple = typename soa_detail::TupleToVectors<FieldTuple>::type;
    using PtrsTuple   = typename soa_detail::TupleToVectorPtrs<FieldTuple>::type;

    /// Staging row -- fill fields here, then call push_back().
    /// Mirrors ScalarFieldsBuffer::data.
    Struct row{};

    // ------------------------------------------------------------------
    // Construction
    // ------------------------------------------------------------------
    VectorFieldsBuffer() : ptrs_(make_ptrs(std::make_index_sequence<kNFields>{})) {}

    explicit VectorFieldsBuffer(std::size_t reserve_n) : VectorFieldsBuffer() {
        reserve(reserve_n);
    }

    // Copying would leave ptrs_ pointing at the source's arrays_ -- disallow.
    VectorFieldsBuffer(const VectorFieldsBuffer&)            = delete;
    VectorFieldsBuffer& operator=(const VectorFieldsBuffer&) = delete;

    // Move is safe: re-initialise ptrs_ to point at the new arrays_.
    VectorFieldsBuffer(VectorFieldsBuffer&& o) noexcept
        : arrays_(std::move(o.arrays_))
        , ptrs_(make_ptrs(std::make_index_sequence<kNFields>{})) {}

    VectorFieldsBuffer& operator=(VectorFieldsBuffer&& o) noexcept {
        if (this != &o) {
            arrays_ = std::move(o.arrays_);
            ptrs_   = make_ptrs(std::make_index_sequence<kNFields>{});
        }
        return *this;
    }

    // ------------------------------------------------------------------
    // Staging row access
    // ------------------------------------------------------------------

    Struct* operator->() noexcept       { return &row; }
    const Struct* operator->() const noexcept { return &row; }

    // ------------------------------------------------------------------
    // Enable / disable
    // ------------------------------------------------------------------

    /// Activate or deactivate write operations (default: enabled).
    /// When disabled, make_branches, push_back(), commit_and_reset(),
    /// and clear() are all no-ops.  Must be called before make_branches.
    void enable(bool e = true) noexcept { enabled_ = e; }

    [[nodiscard]] bool is_enabled() const noexcept { return enabled_; }
    [[nodiscard]] explicit operator bool() const noexcept { return enabled_; }

    // ------------------------------------------------------------------
    // Write path
    // ------------------------------------------------------------------

    /// Commit the current staging row to storage.  No-op when disabled.
    void push_back() {
        if (enabled_) push_back_direct(row);
    }

    /// Append s directly, bypassing the staging row.  Always active.
    void push_back(const Struct& s) {
        push_back_direct(s);
    }

    /// Commit staging row then zero-initialise it.  No-op when disabled.
    void commit_and_reset() {
        if (enabled_) { push_back_direct(row); reset_row(); }
    }

    /// Zero-initialise the staging row without touching storage.
    void reset_row() { row = Struct{}; }

    // ------------------------------------------------------------------
    // Storage
    // ------------------------------------------------------------------

    /// Reconstruct a struct from stored row i.
    /// Throws std::out_of_range if i >= size().
    Struct get(std::size_t i) const {
        if (i >= size())
            throw std::out_of_range(
                "VectorFieldsBuffer::get: index " + std::to_string(i) +
                " out of range (size=" + std::to_string(size()) + ")");
        return get_impl(i, std::make_index_sequence<kNFields>{});
    }

    /// Number of committed rows.
    std::size_t size() const {
        return std::get<0>(arrays_).size();
    }

    /// Direct access to the i-th column vector (type-safe via index).
    template<std::size_t I>
    auto& column() { return std::get<I>(arrays_); }

    template<std::size_t I>
    const auto& column() const { return std::get<I>(arrays_); }

    /// Pre-allocate all column vectors.
    void reserve(std::size_t n) {
        std::apply([n](auto&... vecs) { (vecs.reserve(n), ...); }, arrays_);
    }

    /// Empty storage and zero-initialise the staging row.
    /// No-op when disabled.
    void clear() {
        if (enabled_) {
            std::apply([](auto&... vecs) { (vecs.clear(), ...); }, arrays_);
            reset_row();
        }
    }

    // ------------------------------------------------------------------
    // ROOT TTree interface
    // ------------------------------------------------------------------

    /// Create one STL-vector branch per field.  No-op when disabled.
    void make_branches(TTree& tree, const std::string& prefix = "") {
        if (!enabled_) return;
        auto names = trg_detail::get_field_names<Struct>();
        std::size_t i = 0;
        std::apply([&](auto&... vecs) {
            ((tree.Branch((prefix + names[i++]).c_str(), &vecs)), ...);
        }, arrays_);
    }

    /// Re-point branch addresses to this buffer's vectors.
    void set_branch_addresses(TTree& tree, const std::string& prefix = "") {
        auto names = trg_detail::get_field_names<Struct>();
        std::size_t i = 0;
        std::apply([&](auto&... ptrs) {
            ([&](auto& ptr) {
                const auto name = prefix + names[i++];
                check_set_address(tree.SetBranchAddress(name.c_str(), &ptr), name);
            }(ptrs), ...);
        }, ptrs_);
    }

    // ------------------------------------------------------------------
    // Utilities
    // ------------------------------------------------------------------

    static std::array<std::string, kNFields> field_names() {
        return trg_detail::get_field_names<Struct>();
    }

    void print_summary(std::ostream& os = std::cout) const { os << *this; }

    friend std::ostream& operator<<(std::ostream& os, const VectorFieldsBuffer& buf) {
        auto names = trg_detail::get_field_names<Struct>();
        os << "VectorFieldsBuffer<" << typeid(Struct).name()
           << ">  rows=" << buf.size()
           << "  fields=" << kNFields
           << "  enabled=" << std::boolalpha << buf.enabled_ << '\n';
        std::size_t i = 0;
        std::apply([&](const auto&... vecs) {
            ([&](const auto& vec) {
                os << "  " << names[i++] << ": [";
                for (std::size_t j = 0; j < vec.size(); ++j) {
                    if (j) os << ", ";
                    os << vec[j];
                }
                os << "]\n";
            }(vecs), ...);
        }, buf.arrays_);
        return os;
    }

private:
    bool        enabled_ = true;
    ArraysTuple arrays_;
    PtrsTuple   ptrs_;

    template<std::size_t... Is>
    PtrsTuple make_ptrs(std::index_sequence<Is...>) {
        return PtrsTuple{ &std::get<Is>(arrays_)... };
    }

    void push_back_direct(const Struct& s) {
        push_impl(s, std::make_index_sequence<kNFields>{});
    }

    template<std::size_t... Is>
    void push_impl(const Struct& s, std::index_sequence<Is...>) {
        (std::get<Is>(arrays_).push_back(boost::pfr::get<Is>(s)), ...);
    }

    template<std::size_t... Is>
    Struct get_impl(std::size_t i, std::index_sequence<Is...>) const {
        Struct s{};
        ((boost::pfr::get<Is>(s) = std::get<Is>(arrays_)[i]), ...);
        return s;
    }

    static void check_set_address(Int_t status, const std::string& branch_name) {
        if (status < 0)
            throw std::runtime_error(
                "VectorFieldsBuffer::set_branch_addresses: failed to bind branch \"" +
                branch_name + "\" (ROOT error code " + std::to_string(status) + ")");
    }
};

#endif // VECTOR_FIELDS_BUFFER_HH
