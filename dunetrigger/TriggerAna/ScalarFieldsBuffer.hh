#ifndef SCALAR_FIELDS_BUFFER_HH
#define SCALAR_FIELDS_BUFFER_HH
/**
 * @file ScalarFieldsBuffer.hh
 * @brief Saves one struct per ROOT event as scalar branches (no std::vector).
 *
 * Intended for per-event quantities: run number, event ID, trigger flags,
 * global kinematics, etc.  For per-object collections use VectorFieldsBuffer.hh.
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
#include <string>
#include <string_view>
#include <type_traits>
#include <typeinfo>

/**
 * @brief Holds one instance of @p Struct and registers each field as a
 *        scalar branch on a TTree.
 *
 * The struct is the branch buffer itself — ROOT reads/writes directly into
 * its fields.  Intended for one Fill() per event.
 *
 * @tparam Struct A default-constructible, copy-assignable aggregate whose
 *                fields become scalar ROOT TTree branches.
 */
template<trg_concepts::PfrAggregate Struct>
    requires std::is_default_constructible_v<Struct>
          && std::is_copy_assignable_v<Struct>
class ScalarFieldsBuffer {
public:
    /// @brief Number of fields in @p Struct, determined at compile time.
    static constexpr std::size_t kNFields = boost::pfr::tuple_size_v<Struct>;

    /// @brief The live struct — ROOT branches point directly into its fields.
    Struct data{};

    /// @brief Access fields of the live struct directly.
    Struct* operator->() noexcept { return &data; }
    /// @copydoc operator->()
    const Struct* operator->() const noexcept { return &data; }

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

    /// @copydoc is_enabled()
    [[nodiscard]] explicit operator bool() const noexcept { return enabled_; }

    // ------------------------------------------------------------------
    // Reset
    // ------------------------------------------------------------------

    /// @brief Zero-initialise all fields.
    void reset() { data = Struct{}; }

    // ------------------------------------------------------------------
    // ROOT TTree interface — write
    // ------------------------------------------------------------------

    /// @brief Register one scalar branch per field on @p tree.
    /// @param tree   The TTree to attach branches to.
    /// @param prefix Optional prefix prepended to each branch name.
    /// @note No-op when disabled.  Branch type is auto-deduced from the field type.
    void make_branches(TTree& tree, const std::string& prefix = "") {
        if (!enabled_) return;
        constexpr auto names = trg_detail::get_field_names_sv<Struct>();
        trg_detail::for_fields<kNFields>([&]<std::size_t I>() {
            tree.Branch((prefix + names[I]).c_str(), &boost::pfr::get<I>(data));
        });
    }

    // ------------------------------------------------------------------
    // ROOT TTree interface — read
    // ------------------------------------------------------------------

    /// @brief Point each branch address at the corresponding field of #data.
    /// @param tree   The TTree to read branch addresses from.
    /// @param prefix Optional prefix prepended to each branch name.
    void set_branch_addresses(TTree& tree, const std::string& prefix = "") {
        constexpr auto names = trg_detail::get_field_names_sv<Struct>();
        trg_detail::for_fields<kNFields>([&]<std::size_t I>() {
            tree.SetBranchAddress((prefix + names[I]).c_str(), &boost::pfr::get<I>(data));
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

    /// @brief Print a summary of field names and current values to @p os.
    void print_summary(std::ostream& os = std::cout) const { os << *this; }

    /// @brief Stream insertion operator — prints field names and current values.
    friend std::ostream& operator<<(std::ostream& os, const ScalarFieldsBuffer& buf) {
        constexpr auto names = trg_detail::get_field_names_sv<Struct>();
        os << "ScalarFieldsBuffer<" << typeid(Struct).name()
           << ">  fields=" << kNFields
           << "  enabled=" << std::boolalpha << buf.enabled_ << '\n';
        trg_detail::for_fields<kNFields>([&]<std::size_t I>() {
            os << "  " << names[I] << ": " << boost::pfr::get<I>(buf.data) << '\n';
        });
        return os;
    }

private:
    bool enabled_ = true;
};

#endif // SCALAR_FIELDS_BUFFER_HH
