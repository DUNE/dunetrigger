#ifndef VECTOR_FIELDS_BUFFER_HH
#define VECTOR_FIELDS_BUFFER_HH
/**
 * @file VectorFieldsBuffer.hh
 * @brief Structure-of-Arrays buffer auto-reflected from a C++ struct,
 *        with ROOT TTree branch registration and clear support.
 *
 * @par Requirements  C++20, Boost >= 1.80 (Boost.PFR), ROOT >= 6.x
 *
 * @par Field-name reflection
 * Real field names are derived automatically at compile time via
 * `boost::pfr::names_as_array`.  No registration macros required.
 */

#include "FieldNames.hh"

#include <boost/pfr.hpp>
#include <TTree.h>

#include <array>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <typeinfo>
#include <vector>

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
template<trg_concepts::PfrAggregate Struct>
    requires std::is_default_constructible_v<Struct>
          && std::is_copy_constructible_v<Struct>
          && (boost::pfr::tuple_size_v<Struct> > 0)
class VectorFieldsBuffer {
public:
    /// @brief Number of fields in @p Struct, determined at compile time.
    static constexpr std::size_t kNFields = boost::pfr::tuple_size_v<Struct>;

    /// @brief Tuple of field types mirroring @p Struct.
    using FieldTuple  = decltype(boost::pfr::structure_to_tuple(std::declval<Struct>()));
    /// @brief Tuple of `std::vector` column storage types.
    using ArraysTuple = typename trg_detail::Vectorize<FieldTuple>::arrays;
    /// @brief Tuple of `std::vector*` pointer types used by ROOT's SetBranchAddress.
    using PtrsTuple   = typename trg_detail::Vectorize<FieldTuple>::ptrs;

    /// @brief Staging row — fill fields here, then call push_back().
    Struct row{};

    // ------------------------------------------------------------------
    // Construction
    // ------------------------------------------------------------------

    VectorFieldsBuffer() : ptrs_(rebind_ptrs()) {}

    /// @brief Construct and pre-allocate column vectors.
    /// @param reserve_n Number of rows to reserve in each column vector.
    explicit VectorFieldsBuffer(std::size_t reserve_n) : VectorFieldsBuffer() {
        reserve(reserve_n);
    }

    /// @note Copy-disabled: copying would leave ptrs_ pointing at the source's arrays_.
    VectorFieldsBuffer(const VectorFieldsBuffer&)            = delete;
    /// @note Copy-disabled: copying would leave ptrs_ pointing at the source's arrays_.
    VectorFieldsBuffer& operator=(const VectorFieldsBuffer&) = delete;

    /// @brief Move constructor — re-points ptrs_ at the new arrays_.
    VectorFieldsBuffer(VectorFieldsBuffer&& o) noexcept
        : arrays_(std::move(o.arrays_)), ptrs_(rebind_ptrs()) {}

    /// @brief Move assignment — re-points ptrs_ at the new arrays_.
    VectorFieldsBuffer& operator=(VectorFieldsBuffer&& o) noexcept {
        if (this != &o) {
            arrays_ = std::move(o.arrays_);
            ptrs_   = rebind_ptrs();
        }
        return *this;
    }

    // ------------------------------------------------------------------
    // Staging row access
    // ------------------------------------------------------------------

    /// @brief Access fields of the staging row directly.
    Struct* operator->() noexcept             { return &row; }
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

    /// @copydoc is_enabled()
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
    void push_back(const Struct& s) { push_back_direct(s); }

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
    /// @throws std::out_of_range if @p i >= size().
    Struct get(std::size_t i) const {
        if (i >= size())
            throw std::out_of_range(
                "VectorFieldsBuffer::get: index " + std::to_string(i) +
                " out of range (size=" + std::to_string(size()) + ")");
        Struct s{};
        trg_detail::for_fields<kNFields>([&]<std::size_t I>() {
            boost::pfr::get<I>(s) = std::get<I>(arrays_)[i];
        });
        return s;
    }

    /// @return Number of committed rows.
    [[nodiscard]] std::size_t size() const { return std::get<0>(arrays_).size(); }

    /// @brief Direct access to the I-th column vector (type-safe via index).
    template<std::size_t I>
    auto& column() { return std::get<I>(arrays_); }

    /// @copydoc column()
    template<std::size_t I>
    const auto& column() const { return std::get<I>(arrays_); }

    /// @brief Pre-allocate all column vectors.
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
    /// @note No-op when disabled.
    void make_branches(TTree& tree, const std::string& prefix = "") {
        if (!enabled_) return;
        constexpr auto names = trg_detail::get_field_names_sv<Struct>();
        trg_detail::for_fields<kNFields>([&]<std::size_t I>() {
            tree.Branch((prefix + names[I]).c_str(), &std::get<I>(arrays_));
        });
    }

    /// @brief Re-point branch addresses to this buffer's vectors.
    /// @throws std::runtime_error if a branch cannot be bound.
    void set_branch_addresses(TTree& tree, const std::string& prefix = "") {
        constexpr auto names = trg_detail::get_field_names_sv<Struct>();
        trg_detail::for_fields<kNFields>([&]<std::size_t I>() {
            const auto name = prefix + names[I];
            check_set_address(tree.SetBranchAddress(name.c_str(), &std::get<I>(ptrs_)), name);
        });
    }

    // ------------------------------------------------------------------
    // Utilities
    // ------------------------------------------------------------------

    /// @brief Returns the compile-time array of field name string_views.
    [[nodiscard]] static constexpr std::array<std::string_view, kNFields>
    field_names() noexcept {
        return trg_detail::get_field_names_sv<Struct>();
    }

    /// @brief Print a summary of field names and stored values to @p os.
    void print_summary(std::ostream& os = std::cout) const { os << *this; }

    /// @brief Stream insertion operator — prints field names and stored column contents.
    friend std::ostream& operator<<(std::ostream& os, const VectorFieldsBuffer& buf) {
        constexpr auto names = trg_detail::get_field_names_sv<Struct>();
        os << "VectorFieldsBuffer<" << typeid(Struct).name()
           << ">  rows=" << buf.size()
           << "  fields=" << kNFields
           << "  enabled=" << std::boolalpha << buf.enabled_ << '\n';
        trg_detail::for_fields<kNFields>([&]<std::size_t I>() {
            const auto& vec = std::get<I>(buf.arrays_);
            os << "  " << names[I] << ": [";
            for (std::size_t j = 0; j < vec.size(); ++j) {
                if (j) os << ", ";
                os << vec[j];
            }
            os << "]\n";
        });
        return os;
    }

private:
    bool        enabled_ = true;
    ArraysTuple arrays_;
    PtrsTuple   ptrs_;

    // Re-points each element of ptrs_ at the corresponding vector in arrays_.
    // Called after construction and after any move.
    PtrsTuple rebind_ptrs() noexcept {
        PtrsTuple p;
        trg_detail::for_fields<kNFields>([&]<std::size_t I>() {
            std::get<I>(p) = &std::get<I>(arrays_);
        });
        return p;
    }

    void push_back_direct(const Struct& s) {
        trg_detail::for_fields<kNFields>([&]<std::size_t I>() {
            std::get<I>(arrays_).push_back(boost::pfr::get<I>(s));
        });
    }

    static void check_set_address(Int_t status, const std::string& branch_name) {
        if (status < 0)
            throw std::runtime_error(
                "VectorFieldsBuffer::set_branch_addresses: failed to bind branch \"" +
                branch_name + "\" (ROOT error code " + std::to_string(status) + ")");
    }
};

#endif // VECTOR_FIELDS_BUFFER_HH
