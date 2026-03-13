#ifndef VECTOR_FIELDS_BUFFER_HH
#define VECTOR_FIELDS_BUFFER_HH
/**
 * @file VectorFieldsBuffer.hh
 * @brief Automatic Structure-of-Arrays buffer from a C++ struct,
 *        with ROOT TTree branch registration and clear support.
 *
 * @par Requirements
 * - C++17 or later 
 * - Boost >= 1.75 (Boost.PFR, header-only)
 *
 * @par Field-name reflection
 * In C++20 mode (Boost >= 1.80) real struct field names are used for
 * branch names (e.g. `trk_px`).  In C++17 mode use
 * `REGISTER_FIELD_NAMES(StructType, field1, ...)` at namespace scope.
 *
 */

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

/** @internal Internal implementation helpers — not part of the public API. */
namespace soa_detail {

/// @brief Maps `tuple<T0, T1, ...>` to `tuple<vector<T0>, vector<T1>, ...>`.
template<typename Tuple>
struct TupleToVectors;

template<typename... Ts>
struct TupleToVectors<std::tuple<Ts...>> {
    using type = std::tuple<std::vector<Ts>...>;
};

/// @brief Maps `tuple<T0, T1, ...>` to `tuple<vector<T0>*, vector<T1>*, ...>`.
/// @note Used to provide the `T**` pointers that ROOT's SetBranchAddress requires.
template<typename Tuple>
struct TupleToVectorPtrs;

template<typename... Ts>
struct TupleToVectorPtrs<std::tuple<Ts...>> {
    using type = std::tuple<std::vector<Ts>*...>;
};

} // namespace soa_detail

/**
 * @brief Wraps one `std::vector<T>` per field of @p Struct, reflecting the
 *        layout automatically through Boost.PFR.
 *
 * @par Write path
 * Fill the public staging row, then commit it:
 * @code
 * buf->field = value;    // or buf.row.field = value
 * buf.push_back();       // commit row to SoA storage
 * tree.Fill();
 * buf.clear();           // empty storage and reset row
 * @endcode
 *
 * @par Read path
 * Attach to an existing tree, then iterate:
 * @code
 * buf.set_branch_addresses(tree, "prefix_");
 * for (Long64_t e = 0; e < tree.GetEntries(); ++e) {
 *     tree.GetEntry(e);
 *     buf.get(i);        // reconstruct one AoS element
 * }
 * @endcode
 *
 * @tparam Struct A default-constructible, copy-constructible aggregate whose
 *                fields are stored as parallel `std::vector` columns.
 */
template<typename Struct>
class VectorFieldsBuffer {
    static_assert(std::is_default_constructible_v<Struct>,
                  "VectorFieldsBuffer requires a default-constructible struct");
    static_assert(std::is_copy_constructible_v<Struct>,
                  "VectorFieldsBuffer requires a copy-constructible struct");
    static_assert(boost::pfr::tuple_size_v<Struct> > 0,
                  "VectorFieldsBuffer does not support empty structs");

public:
    /// @brief Number of fields in @p Struct, determined at compile time.
    static constexpr std::size_t kNFields = boost::pfr::tuple_size_v<Struct>;

    /// @brief Tuple of field types mirroring @p Struct.
    using FieldTuple  = decltype(boost::pfr::structure_to_tuple(std::declval<Struct>()));
    /// @brief Tuple of `std::vector` column storage types.
    using ArraysTuple = typename soa_detail::TupleToVectors<FieldTuple>::type;
    /// @brief Tuple of `std::vector*` pointer types used by ROOT's SetBranchAddress.
    using PtrsTuple   = typename soa_detail::TupleToVectorPtrs<FieldTuple>::type;

    /// @brief Staging row — fill fields here, then call push_back().
    /// @note Mirrors ScalarFieldsBuffer::data.
    Struct row{};

    // ------------------------------------------------------------------
    // Construction
    // ------------------------------------------------------------------

    /// @brief Default constructor.
    VectorFieldsBuffer() : ptrs_(make_ptrs(std::make_index_sequence<kNFields>{})) {}

    /// @brief Construct and pre-allocate column vectors.
    /// @param reserve_n Number of rows to reserve in each column vector.
    explicit VectorFieldsBuffer(std::size_t reserve_n) : VectorFieldsBuffer() {
        reserve(reserve_n);
    }

    /// @note Copy-disabled: copying would leave ptrs_ pointing at the source's arrays_.
    VectorFieldsBuffer(const VectorFieldsBuffer&)            = delete;
    /// @note Copy-disabled: copying would leave ptrs_ pointing at the source's arrays_.
    VectorFieldsBuffer& operator=(const VectorFieldsBuffer&) = delete;

    /// @brief Move constructor — re-initialises ptrs_ to point at the new arrays_.
    VectorFieldsBuffer(VectorFieldsBuffer&& o) noexcept
        : arrays_(std::move(o.arrays_))
        , ptrs_(make_ptrs(std::make_index_sequence<kNFields>{})) {}

    /// @brief Move assignment — re-initialises ptrs_ to point at the new arrays_.
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

    /// @brief Access fields of the staging row directly.
    Struct* operator->() noexcept       { return &row; }
    /// @copydoc operator->()
    const Struct* operator->() const noexcept { return &row; }

    // ------------------------------------------------------------------
    // Enable / disable
    // ------------------------------------------------------------------

    /// @brief Activate or deactivate write operations.
    /// @param e @c true to enable (default), @c false to disable.
    /// @note When disabled, make_branches(), push_back(), commit_and_reset(),
    ///       and clear() are all no-ops.  Must be called before make_branches().
    void enable(bool e = true) noexcept { enabled_ = e; }

    /// @return @c true if write operations are active.
    [[nodiscard]] bool is_enabled() const noexcept { return enabled_; }

    /// @return @c true if write operations are active.
    [[nodiscard]] explicit operator bool() const noexcept { return enabled_; }

    // ------------------------------------------------------------------
    // Write path
    // ------------------------------------------------------------------

    /// @brief Commit the current staging row to storage.
    /// @note No-op when disabled.
    void push_back() {
        if (enabled_) push_back_direct(row);
    }

    /// @brief Append @p s directly, bypassing the staging row.
    /// @param s Row to append.  Always active regardless of the enabled flag.
    void push_back(const Struct& s) {
        push_back_direct(s);
    }

    /// @brief Commit the staging row then zero-initialise it.
    /// @note No-op when disabled.
    void commit_and_reset() {
        if (enabled_) { push_back_direct(row); reset_row(); }
    }

    /// @brief Zero-initialise the staging row without touching storage.
    void reset_row() { row = Struct{}; }

    // ------------------------------------------------------------------
    // Storage
    // ------------------------------------------------------------------

    /// @brief Reconstruct a struct from stored row @p i.
    /// @param i Zero-based row index.
    /// @return Copy of the stored row at index @p i.
    /// @throws std::out_of_range if @p i >= size().
    Struct get(std::size_t i) const {
        if (i >= size())
            throw std::out_of_range(
                "VectorFieldsBuffer::get: index " + std::to_string(i) +
                " out of range (size=" + std::to_string(size()) + ")");
        return get_impl(i, std::make_index_sequence<kNFields>{});
    }

    /// @return Number of committed rows.
    std::size_t size() const {
        return std::get<0>(arrays_).size();
    }

    /// @brief Direct access to the I-th column vector (type-safe via index).
    /// @tparam I Zero-based column index.
    /// @return Reference to the column's `std::vector`.
    template<std::size_t I>
    auto& column() { return std::get<I>(arrays_); }

    /// @copydoc column()
    template<std::size_t I>
    const auto& column() const { return std::get<I>(arrays_); }

    /// @brief Pre-allocate all column vectors.
    /// @param n Number of rows to reserve.
    void reserve(std::size_t n) {
        std::apply([n](auto&... vecs) { (vecs.reserve(n), ...); }, arrays_);
    }

    /// @brief Empty storage and zero-initialise the staging row.
    /// @note No-op when disabled.
    void clear() {
        if (enabled_) {
            std::apply([](auto&... vecs) { (vecs.clear(), ...); }, arrays_);
            reset_row();
        }
    }

    // ------------------------------------------------------------------
    // ROOT TTree interface
    // ------------------------------------------------------------------

    /// @brief Create one STL-vector branch per field on @p tree.
    /// @param tree   The TTree to attach branches to.
    /// @param prefix Optional prefix prepended to each branch name.
    /// @note No-op when disabled.
    void make_branches(TTree& tree, const std::string& prefix = "") {
        if (!enabled_) return;
        auto names = trg_detail::get_field_names<Struct>();
        std::size_t i = 0;
        std::apply([&](auto&... vecs) {
            ((tree.Branch((prefix + names[i++]).c_str(), &vecs)), ...);
        }, arrays_);
    }

    /// @brief Re-point branch addresses to this buffer's vectors.
    /// @param tree   The TTree to bind branch addresses from.
    /// @param prefix Optional prefix prepended to each branch name.
    /// @throws std::runtime_error if a branch cannot be bound.
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

    /// @brief Returns the array of field names for @p Struct.
    /// @return Array of @c kNFields field-name strings.
    static std::array<std::string, kNFields> field_names() {
        return trg_detail::get_field_names<Struct>();
    }

    /// @brief Print a summary of field names and stored values to @p os.
    /// @param os Output stream (default: @c std::cout).
    void print_summary(std::ostream& os = std::cout) const { os << *this; }

    /// @brief Stream insertion operator — prints field names and stored column contents.
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
