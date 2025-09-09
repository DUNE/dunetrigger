/**
 * @file TriggerActivityMakerChargeDensity.hpp
 *
 * This is part of the DUNE DAQ Application Framework, copyright 2021.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */

#ifndef TRIGGERALGS_CHARGEDENSITY_TRIGGERACTIVITYMAKERCHARGEDENSITY_HPP_
#define TRIGGERALGS_CHARGEDENSITY_TRIGGERACTIVITYMAKERCHARGEDENSITY_HPP_

#include "dunetrigger/triggeralgs/include/triggeralgs/TriggerActivityFactory.hpp"
#include "dunetrigger/triggeralgs/include/triggeralgs/Types.hpp"

#include <optional>
#include <vector>

namespace triggeralgs {
class TriggerActivityMakerChargeDensity : public TriggerActivityMaker {

public:
  void operator()(const TriggerPrimitive &input_tp,
                  std::vector<TriggerActivity> &output_ta);

  void configure(const nlohmann::json &config);

  std::optional<TriggerActivity>
  construct_ta(const std::vector<TriggerPrimitive> &tps) const;

  uint64_t m_primitive_count = 0;

  // Configurable parameters.
  uint32_t m_channel_window = 6;
  timestamp_t m_time_window = 1200;
  uint32_t m_adc_threshold = 1200000;
  std::map<uint32_t, std::vector<TriggerPrimitive>> channel_bins;
  uint64_t curr_time_bin;
};
} // namespace triggeralgs

#endif // TRIGGERALGS_CHARGEDENSITY_TRIGGERACTIVITYMAKERCHARGEDENSITY_HPP_
