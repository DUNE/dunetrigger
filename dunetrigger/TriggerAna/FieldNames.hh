#ifndef TRG_FIELD_NAMES_HH
#define TRG_FIELD_NAMES_HH
// =============================================================================
//  FieldNames.hpp
//  Shared field-name reflection used by SoABuffer and ScalarBuffer.
//
//  Provides:
//    FieldNames<Struct>               -- specialisable trait
//    REGISTER_FIELD_NAMES(Type, ...)  -- macro to register names (C++17)
//    trg_detail::get_field_names<S>() -- runtime std::array of branch names
//
//  C++20: names are derived automatically via boost::pfr::names_as_array.
//  C++17: call REGISTER_FIELD_NAMES at namespace scope after the struct.
// =============================================================================

#include <boost/pfr.hpp>

#include <array>
#include <string>
#include <utility>

#include <boost/preprocessor/variadic/to_seq.hpp>
#include <boost/preprocessor/seq/transform.hpp>
#include <boost/preprocessor/seq/enum.hpp>
#include <boost/preprocessor/stringize.hpp>

// ---------------------------------------------------------------------------
//  FieldNames<Struct>
//  Specialisable trait.  Default produces "field_0", "field_1", ...
// ---------------------------------------------------------------------------
template<typename Struct>
struct FieldNames {
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
//  Boost.PP stringify helpers -- internal, prefixed TRG_ to avoid clashes.
// ---------------------------------------------------------------------------
#define TRG_PP_STRINGIFY_OP(r, _, elem)  BOOST_PP_STRINGIZE(elem)
#define TRG_PP_STRINGIFY_EACH(...)                              \
    BOOST_PP_SEQ_ENUM(                                          \
        BOOST_PP_SEQ_TRANSFORM(                                 \
            TRG_PP_STRINGIFY_OP, _,                             \
            BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__)))

// ---------------------------------------------------------------------------
//  REGISTER_FIELD_NAMES(StructType, field1, field2, ...)
//  Specialises FieldNames<StructType>.  Place at namespace scope after the
//  struct definition.  Enforces that the name count matches the field count.
// ---------------------------------------------------------------------------
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
};

// ---------------------------------------------------------------------------
//  trg_detail::get_field_names<Struct>()
//  Returns a runtime std::array of branch-name strings.
//  C++20: automatic via boost::pfr::names_as_array.
//  C++17: delegates to FieldNames<Struct>::get() (must be registered).
// ---------------------------------------------------------------------------
namespace trg_detail {

#if __cplusplus >= 202002L
// Forward-declare before get_field_names to satisfy all compilers.
template<typename Struct, typename NamesArray, std::size_t... Is>
std::array<std::string, sizeof...(Is)>
get_names_impl(const NamesArray& pfr_names, std::index_sequence<Is...>);
#endif

template<typename Struct>
std::array<std::string, boost::pfr::tuple_size_v<Struct>>
get_field_names() {
#if __cplusplus >= 202002L
    constexpr auto pfr_names = boost::pfr::names_as_array<Struct>();
    return get_names_impl<Struct>(pfr_names,
        std::make_index_sequence<boost::pfr::tuple_size_v<Struct>>{});
#else
    static_assert(FieldNames<Struct>::registered,
        "C++17 mode: field names not registered for this struct. "
        "Use REGISTER_FIELD_NAMES(StructType, field1, ...) "
        "at namespace scope, or compile with -std=c++20.");
    return FieldNames<Struct>::get();
#endif
}

#if __cplusplus >= 202002L
template<typename Struct, typename NamesArray, std::size_t... Is>
std::array<std::string, sizeof...(Is)>
get_names_impl(const NamesArray& pfr_names, std::index_sequence<Is...>) {
    return { std::string(pfr_names[Is])... };
}
#endif

} // namespace trg_detail

#endif // TRG_FIELD_NAMES_HH
