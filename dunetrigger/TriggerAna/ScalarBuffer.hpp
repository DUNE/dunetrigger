#ifndef SCALAR_BUFFER_HPP
#define SCALAR_BUFFER_HPP
// =============================================================================
//  ScalarBuffer.hpp
//  Saves one POD struct per ROOT event as scalar branches (no std::vector).
//
//  Intended for per-event quantities: run number, event ID, trigger flags,
//  global kinematics, etc.  For per-object collections use SoABuffer.hpp.
//
//  Requirements:
//    - C++17 or later  (GCC >= 12 supported)
//    - Boost >= 1.75   (Boost.PFR + Boost.Preprocessor, header-only)
//    - ROOT >= 6.x     (for TTree / TBranch)
//
//  Field-name reflection:
//    C++20 : real field names via boost::pfr::names_as_array (automatic).
//    C++17 : use REGISTER_FIELD_NAMES(StructType, field1, field2, ...)
//            at namespace scope (provided by FieldNames.hpp).
//
//  Compile (C++17):
//    g++ -std=c++17 main.cpp $(root-config --cflags --libs) -I/path/to/boost -o demo
//  Compile (C++20, real field names):
//    g++ -std=c++20 main.cpp $(root-config --cflags --libs) -I/path/to/boost -o demo
// =============================================================================

#include "FieldNames.hpp"

#include <boost/pfr.hpp>
#include <TTree.h>

#include <array>
#include <iostream>
#include <string>
#include <type_traits>
#include <utility>


// ---------------------------------------------------------------------------
//  ScalarBuffer<Struct>
//  ---
//  Holds one instance of Struct and registers each field as a scalar branch
//  on a TTree.  The struct is the branch buffer itself -- ROOT reads/writes
//  directly into its fields.  One Fill() per event.
//
//  Methods:
//    make_branches(TTree&, prefix)       -- register scalar branches (write)
//    set_branch_addresses(TTree&, prefix)-- attach to existing branches (read)
//    data                               -- public Struct instance (read/write)
//    reset()                            -- zero-initialise data
//    print_summary(os)                  -- list field names and addresses
// ---------------------------------------------------------------------------
template<typename Struct>
class ScalarBuffer {
    static_assert(std::is_trivially_copyable_v<Struct>,
                  "ScalarBuffer requires a POD / trivially-copyable struct");

public:
    static constexpr std::size_t kNFields = boost::pfr::tuple_size_v<Struct>;

    /// The live struct -- ROOT branches point directly into its fields.
    Struct data{};

    Struct* operator->() noexcept { return &data; }
    const Struct* operator->() const noexcept { return &data; }

    // ------------------------------------------------------------------
    // Construction
    // ------------------------------------------------------------------
    ScalarBuffer() = default;

    // ------------------------------------------------------------------
    // Reset
    // ------------------------------------------------------------------

    /// Zero-initialise all fields.
    void reset() { data = Struct{}; }

    // ------------------------------------------------------------------
    // ROOT TTree interface -- write
    // ------------------------------------------------------------------

    /// Register one scalar branch per field.
    /// Branch leaf descriptor is derived automatically from the field type.
    /// prefix + field_name is used as the branch name.
    void make_branches(TTree& tree, const std::string& prefix = "") {
        make_branches_impl(tree, prefix, std::make_index_sequence<kNFields>{});
    }

    // ------------------------------------------------------------------
    // ROOT TTree interface -- read
    // ------------------------------------------------------------------

    /// Point each branch address at the corresponding field of data.
    void set_branch_addresses(TTree& tree, const std::string& prefix = "") {
        set_addresses_impl(tree, prefix, std::make_index_sequence<kNFields>{});
    }

    // ------------------------------------------------------------------
    // Utilities
    // ------------------------------------------------------------------

    static std::array<std::string, kNFields> field_names() {
        return trg_detail::get_field_names<Struct>();
    }

    void print_summary(std::ostream& os = std::cout) const {
        auto names = trg_detail::get_field_names<Struct>();
        os << "ScalarBuffer<" << typeid(Struct).name()
           << ">  fields=" << kNFields << '\n';
        print_impl(os, names, std::make_index_sequence<kNFields>{});
    }

private:
    // ---- make_branches ---------------------------------------------------
    template<std::size_t... Is>
    void make_branches_impl(TTree& tree,
                            const std::string& prefix,
                            std::index_sequence<Is...>) {
        auto names = trg_detail::get_field_names<Struct>();
        ((tree.Branch((prefix + names[Is]).c_str(),
                      &boost::pfr::get<Is>(data))), ...);
    }

    // ---- set_branch_addresses --------------------------------------------
    template<std::size_t... Is>
    void set_addresses_impl(TTree& tree,
                            const std::string& prefix,
                            std::index_sequence<Is...>) {
        auto names = trg_detail::get_field_names<Struct>();
        (tree.SetBranchAddress(
            (prefix + names[Is]).c_str(),
            &boost::pfr::get<Is>(data)
        ), ...);
    }

    // ---- print -----------------------------------------------------------
    template<std::size_t... Is>
    void print_impl(std::ostream& os,
                    const std::array<std::string, kNFields>& names,
                    std::index_sequence<Is...>) const {
        ((os << "  [" << Is << "] " << names[Is]
             << "  addr=" << static_cast<const void*>(&boost::pfr::get<Is>(data))
             << '\n'), ...);
    }
};

#endif // SCALAR_BUFFER_HPP
