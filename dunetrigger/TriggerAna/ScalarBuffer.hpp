#ifndef SCALAR_BUFFER_HPP
#define SCALAR_BUFFER_HPP
// =============================================================================
//  ScalarBuffer.hpp
//  Saves one POD struct per ROOT event as scalar branches (no std::vector).
//
//  Intended for per-event quantities: run number, event ID, trigger flags,
//  global kinematics, etc.  For per-object collections use SoABuffer.hpp.
//
//  Requirements: C++20, Boost >= 1.80 (Boost.PFR), ROOT >= 6.x
//
//  Unlike ScalarFieldsBuffer, this class uses explicit ROOT leaf descriptors
//  (e.g. "run/I", "energy/F") via the three-argument TTree::Branch form,
//  which is required for some older ROOT workflows.
// =============================================================================

#include "FieldNames.hh"

#include <boost/pfr.hpp>
#include <TTree.h>

#include <array>
#include <iostream>
#include <string>
#include <string_view>
#include <type_traits>
#include <typeinfo>

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
//  on a TTree.  The struct is the branch buffer itself — ROOT reads/writes
//  directly into its fields.  One Fill() per event.
//
//  Methods:
//    make_branches(TTree&, prefix)        — register scalar branches (write)
//    set_branch_addresses(TTree&, prefix) — attach to existing branches (read)
//    data                                 — public Struct instance (read/write)
//    reset()                              — zero-initialise data
//    field_names()                        — constexpr array of string_views
//    print_summary(os)                    — list field names and addresses
// ---------------------------------------------------------------------------
template<trg_concepts::PfrAggregate Struct>
    requires std::is_trivially_copyable_v<Struct>
class ScalarBuffer {
public:
    static constexpr std::size_t kNFields = boost::pfr::tuple_size_v<Struct>;

    /// The live struct — ROOT branches point directly into its fields.
    Struct data{};

    Struct* operator->() noexcept       { return &data; }
    const Struct* operator->() const noexcept { return &data; }

    ScalarBuffer() = default;

    /// @brief Zero-initialise all fields.
    void reset() { data = Struct{}; }

    // ------------------------------------------------------------------
    // ROOT TTree interface — write
    // ------------------------------------------------------------------

    /// @brief Register one scalar branch per field using explicit leaf descriptors.
    /// @param tree   The TTree to attach branches to.
    /// @param prefix Optional prefix prepended to each branch name.
    void make_branches(TTree& tree, const std::string& prefix = "") {
        constexpr auto names = trg_detail::get_field_names_sv<Struct>();
        trg_detail::for_fields<kNFields>([&]<std::size_t I>() {
            using FieldType = boost::pfr::tuple_element_t<I, Struct>;
            const std::string name      = prefix + std::string(names[I]);
            const std::string leaf_list = name + "/" + scalar_detail::RootLeafType<FieldType>::value;
            tree.Branch(name.c_str(), &boost::pfr::get<I>(data), leaf_list.c_str());
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

    /// @brief Print field names and their addresses in @p data.
    void print_summary(std::ostream& os = std::cout) const {
        constexpr auto names = trg_detail::get_field_names_sv<Struct>();
        os << "ScalarBuffer<" << typeid(Struct).name()
           << ">  fields=" << kNFields << '\n';
        trg_detail::for_fields<kNFields>([&]<std::size_t I>() {
            os << "  [" << I << "] " << names[I]
               << "  addr=" << static_cast<const void*>(&boost::pfr::get<I>(data)) << '\n';
        });
    }
};

#endif // SCALAR_BUFFER_HPP
