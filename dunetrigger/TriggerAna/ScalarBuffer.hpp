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
//    C++17 : use REGISTER_SCALAR_FIELD_NAMES(StructType, field1, field2, ...)
//            at namespace scope -- shares the same Boost.PP machinery as
//            REGISTER_SOA_FIELD_NAMES in SoABuffer.hpp.
//
//  Compile (C++17):
//    g++ -std=c++17 main.cpp $(root-config --cflags --libs) -I/path/to/boost -o demo
//  Compile (C++20, real field names):
//    g++ -std=c++20 main.cpp $(root-config --cflags --libs) -I/path/to/boost -o demo
// =============================================================================

#include <boost/pfr.hpp>
#include <TTree.h>

#include <array>
#include <iostream>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <stdexcept>

#include <boost/preprocessor/variadic/to_seq.hpp>
#include <boost/preprocessor/seq/transform.hpp>
#include <boost/preprocessor/seq/enum.hpp>
#include <boost/preprocessor/stringize.hpp>

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

// Build runtime field name array (shared logic with SoABuffer).
// C++20: real names via pfr::names_as_array.
// C++17: delegated to ScalarFieldNames<Struct> trait (see below).
template<typename Struct>
std::array<std::string, boost::pfr::tuple_size_v<Struct>>
get_field_names();   // defined after ScalarFieldNames

} // namespace scalar_detail

// ---------------------------------------------------------------------------
//  ScalarFieldNames<Struct> -- same pattern as SoaFieldNames in SoABuffer.hpp.
//  Specialise via REGISTER_SCALAR_FIELD_NAMES in C++17 mode.
// ---------------------------------------------------------------------------
template<typename Struct>
struct ScalarFieldNames {
    static constexpr bool registered = false;

    static std::array<std::string, boost::pfr::tuple_size_v<Struct>> get() {
        return get_impl(std::make_index_sequence<boost::pfr::tuple_size_v<Struct>>{});
    }
private:
    template<std::size_t... Is>
    static std::array<std::string, sizeof...(Is)>
    get_impl(std::index_sequence<Is...>) {
        return { ("field_" + std::to_string(Is))... };
    }
};

// ---------------------------------------------------------------------------
//  Boost.PP stringify helper (same technique as in SoABuffer.hpp).
// ---------------------------------------------------------------------------
#define SCALAR_PP_STRINGIFY_OP(r, _, elem) BOOST_PP_STRINGIZE(elem)

#define SCALAR_PP_STRINGIFY_EACH(...)                       \
    BOOST_PP_SEQ_ENUM(                                      \
        BOOST_PP_SEQ_TRANSFORM(                             \
            SCALAR_PP_STRINGIFY_OP, _,                      \
            BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__)))

// ---------------------------------------------------------------------------
//  REGISTER_SCALAR_FIELD_NAMES(StructType, field1, field2, ...)
//  Specialises ScalarFieldNames for StructType.
//  Place at namespace scope, after the struct definition.
// ---------------------------------------------------------------------------
#define REGISTER_SCALAR_FIELD_NAMES(StructType, ...)                            \
template<>                                                                       \
struct ScalarFieldNames<StructType> {                                            \
    static constexpr bool registered = true;                                     \
    static constexpr std::size_t kN = boost::pfr::tuple_size_v<StructType>;     \
    static std::array<std::string, kN> get() {                                  \
        static const char* const names[] = {                                    \
            SCALAR_PP_STRINGIFY_EACH(__VA_ARGS__)                                \
        };                                                                       \
        return get_impl(std::make_index_sequence<kN>{}, names);                  \
    }                                                                            \
private:                                                                         \
    template<std::size_t... Is>                                                  \
    static std::array<std::string, sizeof...(Is)>                               \
    get_impl(std::index_sequence<Is...>, const char* const* n) {                \
        return { std::string(n[Is])... };                                        \
    }                                                                            \
};

// ---------------------------------------------------------------------------
//  scalar_detail::get_field_names  (needs ScalarFieldNames to be visible)
// ---------------------------------------------------------------------------
namespace scalar_detail {

template<typename Struct>
std::array<std::string, boost::pfr::tuple_size_v<Struct>>
get_field_names() {
#if __cplusplus >= 202002L
    constexpr auto pfr_names = boost::pfr::names_as_array<Struct>();
    return get_names_impl<Struct>(pfr_names,
        std::make_index_sequence<boost::pfr::tuple_size_v<Struct>>{});
#else
    static_assert(ScalarFieldNames<Struct>::registered,
        "ScalarBuffer (C++17): field names not registered for this struct. "
        "Use REGISTER_SCALAR_FIELD_NAMES(StructType, field1, field2, ...) "
        "at namespace scope, or compile with -std=c++20.");
    return ScalarFieldNames<Struct>::get();
#endif
}

#if __cplusplus >= 202002L
template<typename Struct, typename NamesArray, std::size_t... Is>
std::array<std::string, sizeof...(Is)>
get_names_impl(const NamesArray& pfr_names, std::index_sequence<Is...>) {
    return { std::string(pfr_names[Is])... };
}
#endif

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
        return scalar_detail::get_field_names<Struct>();
    }

    void print_summary(std::ostream& os = std::cout) const {
        auto names = scalar_detail::get_field_names<Struct>();
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
        auto names = scalar_detail::get_field_names<Struct>();
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
        auto names = scalar_detail::get_field_names<Struct>();
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
