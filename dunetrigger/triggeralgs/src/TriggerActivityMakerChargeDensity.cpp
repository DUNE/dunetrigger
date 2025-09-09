/**
 * @file TriggerActivityMakerChargeDensity.cpp
 *
 * This is part of the DUNE DAQ Application Framework, copyright 2021.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */

#include "dunetrigger/triggeralgs/include/triggeralgs/ChargeDensity/TriggerActivityMakerChargeDensity.hpp"

#include "TRACE/trace.h"
#include <limits>
#define TRACE_NAME "TriggerActivityMakerChargeDensityPlugin"

#include <vector>

using namespace triggeralgs;
using Logging::TLVL_DEBUG_ALL;
using Logging::TLVL_DEBUG_HIGH;
using Logging::TLVL_DEBUG_LOW;
using Logging::TLVL_IMPORTANT;

void TriggerActivityMakerChargeDensity::configure(
    const nlohmann::json &config) {
  // FIXME use some schema here
  if (config.is_object()) {
    if (config.contains("channel_window"))
      m_channel_window = config["channel_window"];
    if (config.contains("time_window"))
      m_time_window = config["time_window"];
    if (config.contains("adc_threshold"))
      m_adc_threshold = config["adc_threshold"];
  } else {
    TLOG_DEBUG(TLVL_IMPORTANT) << "[TAM:ChargeDensity] The DEFAULT values for "
                                  "time and channel window is being used.";
  }
  TLOG_DEBUG(TLVL_IMPORTANT)
      << "[TAM:ChargeDensity] Using Channel window: " << m_channel_window
      << " and Time window " << m_time_window << " ADC Threshold "
      << m_adc_threshold;
}

void TriggerActivityMakerChargeDensity::operator()(
    const TriggerPrimitive &input_tp, std::vector<TriggerActivity> &output_ta) {
  m_primitive_count++;
  uint32_t tp_time_bin = input_tp.time_start / m_time_window;
  uint32_t tp_channel_bin = input_tp.channel / m_channel_window;
  if (tp_time_bin > curr_time_bin) {
    for (auto const &[channel_bin_idx, tps] : channel_bins) {
      std::optional<TriggerActivity> ta = construct_ta(tps);
      if (ta)
        output_ta.push_back(ta.value());
    }
    channel_bins.clear();
    curr_time_bin = tp_time_bin;
  }
  channel_bins[tp_channel_bin].push_back(input_tp);
  return;
}

std::optional<TriggerActivity> TriggerActivityMakerChargeDensity::construct_ta(
    const std::vector<TriggerPrimitive> &tps) const {
  TLOG_DEBUG(TLVL_DEBUG_LOW)
      << "[TAM:ChargeDensity] I am constructing a trigger activity!";

  std::optional<TriggerActivity> ta = std::nullopt;
  uint64_t adc_integral = 0UL;
  for (auto const &tp : tps) {
    adc_integral += tp.adc_integral;
  }
  if (adc_integral > m_adc_threshold) {
    ta = TriggerActivity();
    ta->adc_integral = adc_integral;
    ta->time_start = tps[0].time_start;
    ta->time_end = tps[tps.size() - 1].time_start;
    ta->detid = tps[0].detid;
    ta->type = TriggerActivity::Type::kTPC;
    ta->algorithm =
        TriggerActivity::Algorithm::kADCSimpleWindow; // FIXME: obviously!
    ta->time_activity =
        dunedaq::trgdataformats::INVALID_TIMESTAMP; // not used for this TAM

    ta->channel_start = std::numeric_limits<channel_t>::max();
    ta->channel_end = std::numeric_limits<channel_t>::min();
    for (auto const &tp : tps) {
      ta->channel_start = std::min(tp.channel, ta->channel_start);
      ta->channel_end = std::max(tp.channel, ta->channel_end);
    }
    auto peak_tp = *std::max_element(
        tps.begin(), tps.end(),
        [](const TriggerPrimitive &tp1, const TriggerPrimitive &tp2) {
          return tp1.adc_peak < tp2.adc_peak;
        });
    ta->channel_peak = peak_tp.channel;
    ta->time_peak = peak_tp.time_peak;
    ta->adc_peak = peak_tp.adc_peak;
    ta->inputs = tps;
  }
  return ta;
}

// Register algo in TA Factory
REGISTER_TRIGGER_ACTIVITY_MAKER(TRACE_NAME, TriggerActivityMakerChargeDensity)
