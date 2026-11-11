/**
 * @file TriggerActivityROPMakerPassthrough.cpp
 *
 * This is part of the DUNE DAQ Application Framework, copyright 2020.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */

#include "dunetrigger/experimental/tpctriggeralgs/include/tpctriggeralgs/Passthrough/TriggerActivityROPMakerPassthrough.hpp"

#include "dunetrigger/experimental/tpctriggeralgs/include/tpctriggeralgs/ConfigUtils.hpp"

#include "TRACE/trace.h"
#define TRACE_NAME "TriggerActivityROPMakerPassthroughPlugin"

#include <stdexcept>
#include <vector>

namespace tpctriggeralgs {

using Logging::TLVL_DEBUG_MEDIUM;
using Logging::TLVL_IMPORTANT;

void
TriggerActivityROPMakerPassthrough::operator()(const TriggerPrimitive& input_tp, std::vector<TriggerActivityROP>& output_ta)
{
  if ((m_primitive_count++) % m_prescale == 0) {

    TLOG_DEBUG(TLVL_DEBUG_MEDIUM) << "[TAM:Pr] Emitting prescaled TriggerActivityROP " << (m_primitive_count - 1);
    std::vector<TriggerPrimitive> tp_list;
    tp_list.push_back(input_tp);

    TriggerActivityROP ta;
    ta.time_start = input_tp.time_start;
    ta.time_end = input_tp.time_start + input_tp.time_over_threshold;
    ta.time_peak = input_tp.time_peak;
    ta.time_activity = 0;
    ta.channel_start = input_tp.channel;
    ta.channel_end = input_tp.channel;
    ta.channel_peak = input_tp.channel;
    ta.adc_integral = input_tp.adc_integral;
    ta.adc_peak = input_tp.adc_peak;
    ta.detid = input_tp.detid;
    ta.type = static_cast<TriggerActivityROPData::type_t>(TriggerActivityROPData::Type::kTPC);
    ta.algorithm = static_cast<TriggerActivityROPData::algorithm_t>(TriggerActivityROPData::Algorithm::kPassthrough);

    ta.inputs = tp_list;

    output_ta.push_back(ta);
  }
}

void
TriggerActivityROPMakerPassthrough::configure(const nlohmann::json& config)
{
  check_known_keys(config, { "prescale" }, TRACE_NAME);
  if (config.is_object() && config.contains("prescale")) {
    m_prescale = config["prescale"];
  }
  // m_prescale is used as a modulus
  if (m_prescale == 0)
    throw std::invalid_argument(TRACE_NAME ": prescale must be > 0");
  TLOG_DEBUG(TLVL_IMPORTANT) << "[TAM:Pr] Using activity prescale " << m_prescale;
}

// Register algo in TA Factory
REGISTER_TRIGGER_ACTIVITY_ROP_MAKER(TRACE_NAME, TriggerActivityROPMakerPassthrough)

} // namespace tpctriggeralgs
