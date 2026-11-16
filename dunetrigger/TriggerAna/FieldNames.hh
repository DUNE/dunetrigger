#ifndef TRG_FIELD_NAMES_HH
#define TRG_FIELD_NAMES_HH
// =============================================================================
//  FieldNames.hh
//  Shared field-name reflection used by VectorFieldsBuffer and ScalarFieldsBuffer.
//
//  Provides:
//    FieldNames<Struct>               -- specialisable trait
//    REGISTER_FIELD_NAMES(Type, ...)  -- macro to register names (C++17)
//    trg_detail::get_field_names<S>() -- runtime std::array of branch names
//
//  C++20: names are derived automatically via boost::pfr::names_as_array.
//          Requires Boost >= 1.84 (boost::pfr::names_as_array); checked below.
//  C++17: call REGISTER_FIELD_NAMES inside the struct's namespace (dunetrigger)
//          after the struct.  The macro checks at compile time that the names
//          match the struct members in declaration order and in number, so the
//          registered names equal what Boost.PFR yields in C++20.
// =============================================================================

#include <boost/pfr.hpp>
#include <boost/version.hpp>

#include <array>
#include <cstddef>
#include <initializer_list>
#include <string>
#include <utility>

#include <boost/preprocessor/variadic/to_seq.hpp>
#include <boost/preprocessor/seq/transform.hpp>
#include <boost/preprocessor/seq/enum.hpp>
#include <boost/preprocessor/stringize.hpp>

// ---------------------------------------------------------------------------
//  FieldNames<Struct>
//  Specialisable trait.  The primary template only marks the struct as
//  unregistered; REGISTER_FIELD_NAMES provides the specialisation with get().
// ---------------------------------------------------------------------------
namespace dunetrigger {
template<typename Struct>
struct FieldNames {
    static constexpr bool registered = false;
};
} // namespace dunetrigger

// ---------------------------------------------------------------------------
//  Helpers for the REGISTER_FIELD_NAMES declaration-order check.
// ---------------------------------------------------------------------------
namespace trg_detail {
constexpr bool strictly_increasing(std::initializer_list<std::size_t> o) {
    const std::size_t* p = o.begin();
    for (std::size_t i = 1; i < o.size(); ++i)
        if (p[i] <= p[i-1]) return false;
    return true;
}
constexpr std::size_t count(std::initializer_list<std::size_t> o) { return o.size(); }
} // namespace trg_detail

// ---------------------------------------------------------------------------
//  Boost.PP stringify helpers -- internal, prefixed TRG_ to avoid clashes.
// ---------------------------------------------------------------------------
#define TRG_PP_STRINGIFY_OP(r, _, elem)  BOOST_PP_STRINGIZE(elem)
#define TRG_PP_STRINGIFY_EACH(...)                              \
    BOOST_PP_SEQ_ENUM(                                          \
        BOOST_PP_SEQ_TRANSFORM(                                 \
            TRG_PP_STRINGIFY_OP, _,                             \
            BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__)))

#define TRG_PP_OFFSETOF_OP(r, Type, elem) offsetof(Type, elem)
#define TRG_PP_OFFSETOF_EACH(Type, ...)                                           \
    BOOST_PP_SEQ_ENUM(BOOST_PP_SEQ_TRANSFORM(TRG_PP_OFFSETOF_OP, Type,            \
                      BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__)))

// ---------------------------------------------------------------------------
//  REGISTER_FIELD_NAMES(StructType, field1, field2, ...)
//  Specialises FieldNames<StructType>.  Must be invoked inside namespace dunetrigger,
//  after the struct definition.  Enforces that the name count
//  matches the field count.  Invoke with a trailing semicolon.
//  In C++20 mode field names are derived automatically via Boost.PFR and
//  this macro is ignored (it expands to nothing).
//  C++17 check: offsetof of each name must be strictly increasing and the
//  count must equal the field count.  offsetof on non-standard-layout types
//  (rows with std::string) is only conditionally supported; GCC supports it
//  and merely warns, hence the targeted -Winvalid-offsetof suppression.
// ---------------------------------------------------------------------------
#if __cplusplus >= 202002L
static_assert(BOOST_VERSION >= 108400,
    "FieldNames.hh: the C++20 path needs boost::pfr::names_as_array (Boost >= 1.84). "
    "This stack has an older Boost; build in C++17 mode or upgrade Boost.");
#define REGISTER_FIELD_NAMES(StructType, ...) static_assert(true)
#else
#define REGISTER_FIELD_NAMES(StructType, ...)                                   \
template<>                                                                       \
struct FieldNames<StructType> {                                                  \
    static constexpr bool registered = true;                                     \
    static constexpr std::size_t kN = boost::pfr::tuple_size_v<StructType>;     \
    static std::array<std::string, kN> get() {                                  \
        static const char* const names[] = {                                    \
            TRG_PP_STRINGIFY_EACH(__VA_ARGS__)                                   \
        };                                                                       \
        constexpr std::size_t kProvided = sizeof(names) / sizeof(names[0]);     \
        static_assert(kProvided == kN,                                           \
            "REGISTER_FIELD_NAMES: name count does not match "                  \
            "the number of fields in " #StructType ".");                        \
        return get_impl(std::make_index_sequence<kN>{}, names);                  \
    }                                                                            \
private:                                                                         \
    template<std::size_t... Is>                                                  \
    static std::array<std::string, sizeof...(Is)>                               \
    get_impl(std::index_sequence<Is...>, const char* const* n) {                \
        return { std::string(n[Is])... };                                        \
    }                                                                            \
};                                                                               \
_Pragma("GCC diagnostic push")                                                  \
_Pragma("GCC diagnostic ignored \"-Winvalid-offsetof\"")                        \
static_assert(::trg_detail::strictly_increasing(                                \
    { TRG_PP_OFFSETOF_EACH(StructType, __VA_ARGS__) }),                         \
    "REGISTER_FIELD_NAMES: names are not in declaration order for " #StructType); \
static_assert(boost::pfr::tuple_size_v<StructType> ==                           \
    ::trg_detail::count({ TRG_PP_OFFSETOF_EACH(StructType, __VA_ARGS__) }),     \
    "REGISTER_FIELD_NAMES: name count mismatch for " #StructType);              \
_Pragma("GCC diagnostic pop")
#endif // __cplusplus >= 202002L

// ---------------------------------------------------------------------------
//  trg_detail::get_field_names<Struct>()
//  Returns a runtime std::array of branch-name strings.
//  C++20: automatic via boost::pfr::names_as_array.
//  C++17: delegates to dunetrigger::FieldNames<Struct>::get() (must be registered).
// ---------------------------------------------------------------------------
namespace trg_detail {

#if __cplusplus >= 202002L
template<typename Struct, typename NamesArray, std::size_t... Is>
std::array<std::string, sizeof...(Is)>
get_names_impl(const NamesArray& pfr_names, std::index_sequence<Is...>) {
    return { std::string(pfr_names[Is])... };
}
#endif


template<typename Struct>
std::array<std::string, boost::pfr::tuple_size_v<Struct>>
get_field_names() {
#if __cplusplus >= 202002L
    constexpr auto pfr_names = boost::pfr::names_as_array<Struct>();
    return get_names_impl<Struct>(pfr_names,
        std::make_index_sequence<boost::pfr::tuple_size_v<Struct>>{});
#else
    static_assert(::dunetrigger::FieldNames<Struct>::registered,
        "C++17 mode: field names not registered for this struct. "
        "Use REGISTER_FIELD_NAMES(StructType, field1, ...) "
        "inside namespace dunetrigger after the struct, or compile with -std=c++20.");
    // Guarded so an unregistered struct reports only the static_assert above.
    if constexpr (::dunetrigger::FieldNames<Struct>::registered)
        return ::dunetrigger::FieldNames<Struct>::get();
    else
        return {};
#endif
}

} // namespace trg_detail

#endif // TRG_FIELD_NAMES_HH
