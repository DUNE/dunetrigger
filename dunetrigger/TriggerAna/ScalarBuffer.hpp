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
#include <tuple>
#include <type_traits>
#include <utility>
#include <stdexcept>

// ---------------------------------------------------------------------------
//  ROOT type-leaf string for scalar branches.
//  TTree::Branch for scalars requires the type descriptor, e.g. "run/I".
// ---------------------------------------------------------------------------
namespace scalar_detail {

template<typename T> struct RootLeafType;
template<> struct RootLeafType<float>         { static constexpr const char* value = "F"; };
template<> struct RootLeafType<double>        { static constexpr const char* value = "D"; };
template<> struct RootLeafType<int>           { static constexpr const char* value = "I"; };
template<> struct RootLeafType<unsigned int>  { static constexpr const char* value = "i"; };
template<> struct RootLeafType<long>          { static constexpr const char* value = "L"; };
template<> struct RootLeafType<unsigned long> { static constexpr const char* value = "l"; };
template<> struct RootLeafType<short>         { static constexpr const char* value = "S"; };
template<> struct RootLeafType<bool>          { static constexpr const char* value = "O"; };
template<> struct RootLeafType<char>          { static constexpr const char* value = "B"; };

} // namespace scalar_detail

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
    // Scalar branches use the leaf-list form:
    //   tree.Branch("name", &field, "name/T")
    // where T is the ROOT type code for the field's C++ type.
    template<std::size_t... Is>
    void make_branches_impl(TTree& tree,
                            const std::string& prefix,
                            std::index_sequence<Is...>) {
        auto names = trg_detail::get_field_names<Struct>();
        (make_one_branch<Is>(tree, prefix, names[Is]), ...);
    }

    template<std::size_t I>
    void make_one_branch(TTree& tree,
                         const std::string& prefix,
                         const std::string& name) {
        using FieldType = std::remove_reference_t<
            decltype(boost::pfr::get<I>(std::declval<Struct&>()))>;
        const std::string branch_name = prefix + name;
        const std::string leaf_list   = branch_name + "/"
                                      + scalar_detail::RootLeafType<FieldType>::value;
        tree.Branch(branch_name.c_str(),
                    &boost::pfr::get<I>(data),
                    leaf_list.c_str());
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
