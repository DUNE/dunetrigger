/**
 * @file TriggerActivityData.hpp
 *
 * This header defines the TriggerActivityData struct, which
 * aggregates information about a found set of associate trigger
 * primitives (algorithm used, channels and times involved, etc.). It
 * does *not* include per-trigger-primitive information, which need to
 * be associated with TriggerActivityData in a higher level object.
 *
 * This is part of the DUNE DAQ Application Framework, copyright 2020.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */

#ifndef TRGDATAFORMATS_INCLUDE_TRGDATAFORMATS_EXP_TRIGGERACTIVITYDATA_HPP_
#define TRGDATAFORMATS_INCLUDE_TRGDATAFORMATS_EXP_TRIGGERACTIVITYDATA_HPP_

#include "trgdataformats/Types.hpp"
#include <cstdint>
#include <limits>

namespace dunedaq::trgdataformats::exp {

struct TriggerActivityData
{

  static constexpr version_t s_trigger_activity_version = 2;
  static constexpr channel_t s_invalid_channel = std::numeric_limits<channel_t>::max();

  version_t version = s_trigger_activity_version;
  timestamp_t time_start = TypeDefaults::s_invalid_timestamp;
  timestamp_t time_end = TypeDefaults::s_invalid_timestamp;
  timestamp_t time_peak = TypeDefaults::s_invalid_timestamp;
  timestamp_t time_activity = TypeDefaults::s_invalid_timestamp;
  channel_t channel_start = s_invalid_channel;
  channel_t channel_end = s_invalid_channel;
  channel_t channel_peak = s_invalid_channel;
  uint64_t adc_integral = 0;                 // NOLINT(build/unsigned)
  uint16_t adc_peak = 0;                     // NOLINT(build/unsigned)
  detid_t detid = TypeDefaults::s_invalid_detid;
  uint16_t type = 0;
  uint16_t algorithm = 0;
};
  
} // namespace dunedaq::trgdataformats

// This static_assert is meant to alert the developer to bump the
// version if variables are added or removed
static_assert(
	      dunedaq::trgdataformats::TriggerActivityData::s_trigger_activity_version == 2 &&
	      sizeof(dunedaq::trgdataformats::TriggerActivityData) == 80
	      );


#endif // TRGDATAFORMATS_INCLUDE_TRGDATAFORMATS_EXP_TRIGGERACTIVITYDATA_HPP_
