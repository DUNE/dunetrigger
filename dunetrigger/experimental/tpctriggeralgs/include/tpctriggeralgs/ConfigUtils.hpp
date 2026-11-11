/**
 * @file ConfigUtils.hpp
 *
 * Helpers for validating the JSON configuration passed to tpctriggeralgs makers.
 */

#ifndef tpctriggeralgs_INCLUDE_tpctriggeralgs_CONFIGUTILS_HPP_
#define tpctriggeralgs_INCLUDE_tpctriggeralgs_CONFIGUTILS_HPP_

#include <nlohmann/json.hpp>

#include <algorithm>
#include <initializer_list>
#include <stdexcept>
#include <string>
#include <string_view>

namespace tpctriggeralgs {

// Throw if `config` contains any key not listed in `known`.
// Catches misspelled parameters that would otherwise silently fall back to defaults.
inline void
check_known_keys(const nlohmann::json& config, std::initializer_list<std::string_view> known, std::string_view algorithm)
{
  if (!config.is_object())
    return;
  for (const auto& [key, _] : config.items()) {
    if (std::find(known.begin(), known.end(), key) == known.end()) {
      std::string msg = std::string(algorithm) + ": unknown configuration parameter '" + key + "'. Known parameters:";
      for (auto k : known)
        msg += " " + std::string(k);
      throw std::invalid_argument(msg);
    }
  }
}

// Return config[key] converted to T; throw if `key` is missing.
template<typename T>
T
get_required(const nlohmann::json& config, const std::string& key, std::string_view algorithm)
{
  if (!config.is_object() || !config.contains(key))
    throw std::invalid_argument(std::string(algorithm) + ": missing required configuration parameter '" + key + "'");
  return config.at(key).get<T>();
}

} // namespace tpctriggeralgs

#endif // tpctriggeralgs_INCLUDE_tpctriggeralgs_CONFIGUTILS_HPP_
