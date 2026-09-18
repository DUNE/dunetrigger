#ifndef TRG_FIELD_NAMES_HH
#define TRG_FIELD_NAMES_HH
// =============================================================================
//  FieldNames.hh — C++20 field-reflection utilities shared by all buffers.
//
//  Provides:
//    trg_concepts::PfrAggregate          — concept for reflectable aggregates
//    trg_detail::get_field_names_sv<S>() — constexpr array<string_view, N>
//    trg_detail::for_fields<N>(f)        — call f<I>() for each index I in [0,N)
//    trg_detail::Vectorize<Tuple>        — maps tuple<T...> to SoA storage types
//
//  Requirements: C++20, Boost >= 1.80 (Boost.PFR, header-only)
// =============================================================================

#include <boost/pfr.hpp>

#include <array>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace trg_concepts {

/// Satisfied by aggregate types that Boost.PFR can reflect.
template<typename T>
concept PfrAggregate = std::is_aggregate_v<T>;

} // namespace trg_concepts

namespace trg_detail {

/// @brief Returns field names as a constexpr array of string_views (zero allocation).
template<typename Struct>
[[nodiscard]] constexpr auto get_field_names_sv() noexcept {
    return boost::pfr::names_as_array<Struct>();
}

/// @brief Call f.operator()<I>() for each compile-time index I in [0, N).
/// @par Usage
/// @code
/// trg_detail::for_fields<kNFields>([&]<std::size_t I>() {
///     tree.Branch((prefix + names[I]).c_str(), &boost::pfr::get<I>(data));
/// });
/// @endcode
template<std::size_t N, typename F>
constexpr void for_fields(F&& f) {
    [&]<std::size_t... Is>(std::index_sequence<Is...>) {
        (f.template operator()<Is>(), ...);
    }(std::make_index_sequence<N>{});
}

/// @brief Maps a `tuple<T0,T1,...>` to its SoA storage and pointer types.
/// @par Example
/// @code
/// using FieldTuple  = decltype(boost::pfr::structure_to_tuple(...));
/// using ArraysTuple = typename trg_detail::Vectorize<FieldTuple>::arrays;
/// using PtrsTuple   = typename trg_detail::Vectorize<FieldTuple>::ptrs;
/// @endcode
template<typename Tuple>
struct Vectorize;

template<typename... Ts>
struct Vectorize<std::tuple<Ts...>> {
    using arrays = std::tuple<std::vector<Ts>...>;   ///< parallel vector columns
    using ptrs   = std::tuple<std::vector<Ts>*...>;  ///< raw pointers for ROOT read-back
};

} // namespace trg_detail

#endif // TRG_FIELD_NAMES_HH
