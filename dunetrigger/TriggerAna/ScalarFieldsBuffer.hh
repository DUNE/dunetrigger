#ifndef SCALAR_FIELDS_BUFFER_HH
#define SCALAR_FIELDS_BUFFER_HH
// =============================================================================
//  ScalarFieldsBuffer.hh
//  Saves one struct per ROOT event as scalar branches (no std::vector).
//
//  Intended for per-event quantities: run number, event ID, trigger flags,
//  global kinematics, etc.  For per-object collections use VectorFieldsBuffer.hh.
//
//  Requirements:
//    - C++17 or later
//    - Boost >= 1.75   (Boost.PFR + Boost.Preprocessor, header-only)
//
//  Field-name reflection:
//    C++20 : real field names via boost::pfr::names_as_array (automatic).
//    C++17 : use REGISTER_FIELD_NAMES(StructType, field1, field2, ...)
//            at namespace scope (provided by FieldNames.hpp).
// =============================================================================

#include "FieldNames.hh"

#include <boost/pfr.hpp>
#include <TTree.h>

#include <array>
#include <iostream>
#include <string>
#include <type_traits>
#include <utility>


// ---------------------------------------------------------------------------
//  ScalarFieldsBuffer<Struct>
//  ---
//  Holds one instance of Struct and registers each field as a scalar branch
//  on a TTree.  The struct is the branch buffer itself -- ROOT reads/writes
//  directly into its fields.  One Fill() per event.
//
//  Methods:
//    make_branches(TTree&, prefix)       -- register scalar branches (gated)
//    set_branch_addresses(TTree&, prefix)-- attach to existing branches (read)
//    data                               -- public Struct instance (read/write)
//    reset()                            -- zero-initialise data
//    enable(bool) / is_enabled()        -- enable/disable write operations
//    operator bool()                    -- true when enabled
//    print_summary(os)                  -- list field names and addresses
// ---------------------------------------------------------------------------
template<typename Struct>
class ScalarFieldsBuffer {
    static_assert(std::is_default_constructible_v<Struct>,
                  "ScalarFieldsBuffer requires a default-constructible struct");
    static_assert(std::is_copy_assignable_v<Struct>,
                  "ScalarFieldsBuffer requires a copy-assignable struct");

public:
    static constexpr std::size_t kNFields = boost::pfr::tuple_size_v<Struct>;

    /// The live struct -- ROOT branches point directly into its fields.
    Struct data{};

    Struct* operator->() noexcept { return &data; }
    const Struct* operator->() const noexcept { return &data; }

    // ------------------------------------------------------------------
    // Construction
    // ------------------------------------------------------------------
    ScalarFieldsBuffer() = default;

    // ------------------------------------------------------------------
    // Enable / disable
    // ------------------------------------------------------------------

    /// Activate or deactivate write operations (default: enabled).
    /// When disabled, make_branches is a no-op.
    /// Must be called before make_branches.
    void enable(bool e = true) noexcept { enabled_ = e; }

    [[nodiscard]] bool is_enabled() const noexcept { return enabled_; }
    [[nodiscard]] explicit operator bool() const noexcept { return enabled_; }

    // ------------------------------------------------------------------
    // Reset
    // ------------------------------------------------------------------

    /// Zero-initialise all fields.
    void reset() { data = Struct{}; }

    // ------------------------------------------------------------------
    // ROOT TTree interface -- write
    // ------------------------------------------------------------------

    /// Register one scalar branch per field.  No-op when disabled.
    /// Branch leaf descriptor is derived automatically from the field type.
    /// prefix + field_name is used as the branch name.
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

    /// Point each branch address at the corresponding field of data.
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

    static std::array<std::string, kNFields> field_names() {
        return trg_detail::get_field_names<Struct>();
    }

    void print_summary(std::ostream& os = std::cout) const { os << *this; }

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
