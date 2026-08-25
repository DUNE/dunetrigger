/**
 * @file TriggerActivityMakerProtoDUNEBSMWindow.cpp
 *
 * This is part of the DUNE DAQ Application Framework, copyright 2021.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */

#include "dunetrigger/triggeralgs/include/triggeralgs/ProtoDUNEBSMWindow/TriggerActivityMakerProtoDUNEBSMWindow.hpp"
#include "art/Framework/Principal/Handle.h"

#include "TRACE/trace.h"
#define TRACE_NAME "TriggerActivityMakerProtoDUNEBSMWindowPlugin"

#include <vector>
#include <chrono>

using namespace triggeralgs;
using Logging::TLVL_DEBUG_ALL;
using Logging::TLVL_DEBUG_HIGH;
using Logging::TLVL_DEBUG_LOW;
using Logging::TLVL_IMPORTANT;

void
TriggerActivityMakerProtoDUNEBSMWindow::operator()(const TriggerPrimitive& input_tp, std::vector<TriggerActivity>& output_ta)
{
  
  // The first time operator is called, reset
  // window object.
  if (m_current_window.is_empty()) {
    // Reset window with new TP
    m_current_window.reset(input_tp);
    // Initialise last time an XGBoost prediction was made
    m_last_pred_time = input_tp.time_start;
    // Iterate number of TPs in the window
    m_primitive_count++;
    // first time operator is called set ROP first and last channel
    unsigned int detelement = channelMap->get_element_id_from_offline_channel(input_tp.channel);
    unsigned int plane = channelMap->get_plane_from_offline_channel(input_tp.channel);

    // Are we on the collection plane? Use XGBoost model for collection plane TPs
    // evaluate the sum of the TP charge if on induction planes
    // Induction plane IDs = 0, 1
    // Collection plane ID = 2
    if (plane > 1) m_collection_plane = true;

    if (plane == 3) plane = 2; // PD-HD plane id = 3 is collection plane of APA 1 and 3 in LArSoft

    // Use PlaneInfo object to get the first and last channels on plane
    PlaneInfo plane_info = m_det_plane_map.get_plane_info(m_channel_map_name, detelement, plane);
    m_first_channel = static_cast<channel_t>(plane_info.min_channel);
    m_n_channels_on_plane = static_cast<channel_t>(plane_info.n_channels);

    // If we are in PD-VD and on the collection plane use 'effective' channel mapping for CRPs
    // Effective mapping in PD-VD only for collection plane
    if (m_pdvd_map && m_collection_plane) {
      m_pdvd_eff_channel_mapper = std::make_unique<PDVDEffectiveChannelMap>(plane_info.min_channel, plane_info.n_channels);

      m_first_channel = m_pdvd_eff_channel_mapper->remapCollectionPlaneChannel(m_first_channel);
      m_n_channels_on_plane = m_pdvd_eff_channel_mapper->getNEffectiveChannels();
    }

    return;
  } 

  // If the difference between the current TP's start time and the start of the window
  // is less than the specified window size, add the TP to the window.
  if ((input_tp.time_start - m_current_window.time_start) < m_window_length) {
    TLOG_DEBUG(TLVL_DEBUG_HIGH) << "[TAM:PDBSMW] Window not yet complete, adding the input_tp to the window.";
    m_current_window.add(input_tp);
  }

  // If the addition of the current TP to the window would make it longer
  // than the specified window length, don't add it
  // First, if these are not collection plane TPs, just evaluate the total charge
  // If the total charge on the induction plane crosses a threshold, create a TA
  else if (!m_collection_plane && m_current_window.adc_integral > m_adc_threshold_induction) {
    TLOG_DEBUG(TLVL_DEBUG_LOW) << "[TAM:ADCSW] ADC integral in window is greater than specified threshold.";
    output_ta.push_back(construct_ta());
    TLOG_DEBUG(TLVL_DEBUG_HIGH) << "[TAM:ADCSW] Resetting window with input_tp.";
    m_current_window.reset(input_tp);
  }

  // If the addition of the current TP to the window would make it longer
  // than the specified window length, don't add it
  // Instead go through a series of filters and eventually a XGBoost model to determine whether to create a TA
  else if (m_collection_plane &&
      (m_current_window.time_start - m_last_pred_time) > m_bin_length && // check enough time has passed since last window
      m_current_window.adc_integral > m_adc_threshold_collection && // set a low minimum threshold for the ADC integral sum
      compute_treelite_classification() // XGBoost classifier 
      ) {
    TLOG_DEBUG(TLVL_DEBUG_LOW) << "[TAM:PDBSMW] ADC integral in window is greater than specified threshold.";
    output_ta.push_back(construct_ta());
    TLOG_DEBUG(TLVL_DEBUG_HIGH) << "[TAM:PDBSMW] Resetting window with input_tp.";
    m_current_window.reset(input_tp);
  }
  // If it is not, move the window along.
  else {
    TLOG_DEBUG(TLVL_DEBUG_ALL) << "[TAM:PDBSMW] Window is at required length but adc threshold not met, shifting window along.";
    m_current_window.move(input_tp, m_window_length);
  }
  
  TLOG_DEBUG(TLVL_DEBUG_ALL) << "[TAM:PDBSMW] " << m_current_window;

  m_primitive_count++;

  return;
}

void
TriggerActivityMakerProtoDUNEBSMWindow::configure(const nlohmann::json &config)
{
  //FIXME use some schema here
  if (config.is_object()) {
    if (config.contains("is_pdvd")) m_pdvd_map = config["is_pdvd"];
    if (config.contains("adc_threshold_induction")) m_adc_threshold_induction = config["adc_threshold_induction"];
    if (config.contains("bdt_threshold")) {
      uint64_t int_bdt_threshold = config["bdt_threshold"];
      if (int_bdt_threshold <= 100) m_bdt_threshold = static_cast<float>(int_bdt_threshold * 0.01);
      else if (int_bdt_threshold <= 1000) m_bdt_threshold = static_cast<float>(int_bdt_threshold * 0.001);
      else if (int_bdt_threshold <= 10000) m_bdt_threshold = static_cast<float>(int_bdt_threshold * 0.0001);
      else m_bdt_threshold = static_cast<float>(int_bdt_threshold * 0.01);
    }
  }
  else{
    TLOG_DEBUG(TLVL_IMPORTANT) << "[TAM:PDBSMW] The DEFAULT values of window_length and adc_threshold are being used.";
  }

  // Collection plane adc threshold is fixed based on model training
  // HD and VD models may have different optimum thresholds
  if (m_pdvd_map) {
    m_channel_map_name = "PD2VDTPCChannelMap";
  } else {
    m_channel_map_name = "PD2HDTPCChannelMap";
  }

  
  // Window length in time and number of channel and time bins defined by model training
  // so these are kept constant
  m_bin_length = static_cast<timestamp_t>(m_window_length / m_num_timebins);
 
  channelMap = dunedaq::detchannelmaps::make_map(m_channel_map_name);

  m_compiled_model_interface = std::make_unique<CompiledModelInterface>(nbatch, m_pdvd_map);

  const size_t num_feature = m_compiled_model_interface->GetNumFeatures();

  flat_batched_inputs.resize(num_feature);

  flat_batched_Entries.clear();
  for (size_t i = 0; i < num_feature; ++i) {
    union Entry zero;
    zero.fvalue = 0.0;
    flat_batched_Entries.emplace_back(zero);
  }

}

TriggerActivity
TriggerActivityMakerProtoDUNEBSMWindow::construct_ta() const
{
  TLOG_DEBUG(TLVL_DEBUG_LOW) << "[TAM:PDBSMW] I am constructing a trigger activity!";

  //TriggerPrimitive latest_tp_in_window = m_current_bin.tp_list.back();
  TriggerPrimitive latest_tp_in_window = m_current_window.tp_list.back();
  // The time_peak, time_activity, channel_* and adc_peak fields of this TA are irrelevent
  // for the purpose of this trigger alg.
  TriggerActivity ta;
  ta.time_start = m_current_window.time_start;
  ta.time_end = latest_tp_in_window.time_start + latest_tp_in_window.time_over_threshold;
  ta.time_peak = latest_tp_in_window.time_peak;
  ta.time_activity = latest_tp_in_window.time_peak;
  ta.channel_start = latest_tp_in_window.channel;
  ta.channel_end = latest_tp_in_window.channel;
  ta.channel_peak = latest_tp_in_window.channel;
  ta.adc_integral = m_current_window.adc_integral;
  ta.adc_peak = latest_tp_in_window.adc_peak;
  ta.detid = latest_tp_in_window.detid;
  ta.type = TriggerActivity::Type::kTPC;
  ta.algorithm = TriggerActivity::Algorithm::kUnknown;
  ta.inputs = m_current_window.tp_list;
  return ta;
}

bool TriggerActivityMakerProtoDUNEBSMWindow::compute_treelite_classification() {

  m_last_pred_time = m_current_window.time_start;

  m_current_window.bin_window(
      flat_batched_inputs, 
      m_num_timebins, m_bin_length,
      m_num_chanbins, m_n_channels_on_plane, m_first_channel,
      m_pdvd_eff_channel_mapper, m_pdvd_map
      );

  m_current_window.fill_entry_window(flat_batched_Entries, flat_batched_inputs);

  std::vector<float> result(nbatch, 0.0f);
  
  m_compiled_model_interface->Predict(flat_batched_Entries.data(), result.data());

  return m_compiled_model_interface->Classify(result.data(), m_bdt_threshold);
}

// Register algo in TA Factory
REGISTER_TRIGGER_ACTIVITY_MAKER(TRACE_NAME, TriggerActivityMakerProtoDUNEBSMWindow)
