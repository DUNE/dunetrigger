#ifndef SCALAR_FIELDS_BUFFER_HH
#define SCALAR_FIELDS_BUFFER_HH
/**
 * @file ScalarFieldsBuffer.hh
 * @brief Saves one struct per ROOT event as scalar branches (no std::vector).
 *
 * Intended for per-event quantities: run number, event ID, trigger flags,
 * global kinematics, etc.  For per-object collections use VectorFieldsBuffer.hh.
 *
 * @par Requirements
 * - C++17 or later
 * - Boost >= 1.75 (Boost.PFR + Boost.Preprocessor, header-only)
 *
 * @par Field-name reflection
 * - C++20: real field names via `boost::pfr::names_as_array` (automatic).
 * - C++17: use `REGISTER_FIELD_NAMES(StructType, field1, field2, ...)`
 *   at namespace scope (provided by FieldNames.hpp).
 */

#include "FieldNames.hh"

#include <boost/pfr.hpp>
#include <TTree.h>

#include <array>
#include <iostream>
#include <string>
#include <type_traits>
#include <utility>


/**
 * @brief Holds one instance of @p Struct and registers each field as a
 *        scalar branch on a TTree.
 *
 * The struct is the branch buffer itself -- ROOT reads/writes directly into
 * its fields.  Intended for one Fill() per event.
 *
 * @tparam Struct A default-constructible, copy-assignable aggregate whose
 *                fields become scalar ROOT TTree branches.
 */
template<typename Struct>
class ScalarFieldsBuffer {
    static_assert(std::is_default_constructible_v<Struct>,
                  "ScalarFieldsBuffer requires a default-constructible struct");
    static_assert(std::is_copy_assignable_v<Struct>,
                  "ScalarFieldsBuffer requires a copy-assignable struct");

public:
    /// @brief Number of fields in @p Struct, determined at compile time.
    static constexpr std::size_t kNFields = boost::pfr::tuple_size_v<Struct>;

    /// @brief The live struct — ROOT branches point directly into its fields.
    Struct data{};

    /// @brief Access fields of the live struct directly.
    Struct* operator->() noexcept { return &data; }
    /// @copydoc operator->()
    const Struct* operator->() const noexcept { return &data; }

    // ------------------------------------------------------------------
    // Construction
    // ------------------------------------------------------------------

    /// @brief Default constructor.
    ScalarFieldsBuffer() = default;

    // ------------------------------------------------------------------
    // Enable / disable
    // ------------------------------------------------------------------

    /// @brief Activate or deactivate write operations.
    /// @param e @c true to enable (default), @c false to disable.
    /// @note When disabled, make_branches() is a no-op.
    ///       Must be called before make_branches().
    void enable(bool e = true) noexcept { enabled_ = e; }

    /// @return @c true if write operations are active.
    [[nodiscard]] bool is_enabled() const noexcept { return enabled_; }

    /// @return @c true if write operations are active.
    [[nodiscard]] explicit operator bool() const noexcept { return enabled_; }

    // ------------------------------------------------------------------
    // Reset
    // ------------------------------------------------------------------

    /// @brief Zero-initialise all fields.
    void reset() { data = Struct{}; }

    // ------------------------------------------------------------------
    // ROOT TTree interface -- write
    // ------------------------------------------------------------------

    /// @brief Register one scalar branch per field on @p tree.
    /// @param tree   The TTree to attach branches to.
    /// @param prefix Optional prefix prepended to each branch name.
    /// @note No-op when disabled.  The leaf descriptor is derived
    ///       automatically from the field type.
    void make_branches(TTree& tree, const std::string& prefix = "") {
        if (!enabled_) return;
        auto names = trg_detail::get_field_names<Struct>();
        std::size_t i = 0;
        std::apply([&](auto&... fields) {
            ((tree.Branch((prefix + names[i++]).c_str(), &fields)), ...);
        }, boost::pfr::structure_tie(data));
    }

    // ------------------------------------------------------------------
    // ROOT TTree interface -- read
    // ------------------------------------------------------------------

    /// @brief Point each branch address at the corresponding field of #data.
    /// @param tree   The TTree to read branch addresses from.
    /// @param prefix Optional prefix prepended to each branch name.
    void set_branch_addresses(TTree& tree, const std::string& prefix = "") {
        auto names = trg_detail::get_field_names<Struct>();
        std::size_t i = 0;
        std::apply([&](auto&... fields) {
            ((tree.SetBranchAddress((prefix + names[i++]).c_str(), &fields)), ...);
        }, boost::pfr::structure_tie(data));
    }

    // ------------------------------------------------------------------
    // Utilities
    // ------------------------------------------------------------------

    /// @brief Returns the array of field names for @p Struct.
    /// @return Array of @c kNFields field-name strings.
    static std::array<std::string, kNFields> field_names() {
        return trg_detail::get_field_names<Struct>();
    }

    /// @brief Print a summary of field names and current values to @p os.
    /// @param os Output stream (default: @c std::cout).
    void print_summary(std::ostream& os = std::cout) const { os << *this; }

    /// @brief Stream insertion operator — prints field names and current values.
    friend std::ostream& operator<<(std::ostream& os, const ScalarFieldsBuffer& buf) {
        auto names = trg_detail::get_field_names<Struct>();
        os << "ScalarFieldsBuffer<" << typeid(Struct).name()
           << ">  fields=" << kNFields
           << "  enabled=" << std::boolalpha << buf.enabled_ << '\n';
        std::size_t i = 0;
        std::apply([&](const auto&... fields) {
            ([&](const auto& field) {
                os << "  " << names[i++] << ": " << field << '\n';
            }(fields), ...);
        }, boost::pfr::structure_tie(buf.data));
        return os;
    }

private:
    bool enabled_ = true;
};

#endif // SCALAR_FIELDS_BUFFER_HH
