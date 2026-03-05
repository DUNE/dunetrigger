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
  if(m_current_window.is_empty()){
    m_current_window.reset(input_tp);
    m_last_pred_time = input_tp.time_start;
    m_primitive_count++;
    // first time operator is called set ROP first and last channel
    unsigned int detelement = channelMap->get_element_id_from_offline_channel(input_tp.channel);
    unsigned int plane = channelMap->get_plane_from_offline_channel(input_tp.channel);

    if (plane == 3) plane = 2; // PD-HD plane id = 3 is collection plane of APA 1 and 3 in LArSoft

    PlaneInfo plane_info = m_det_plane_map.get_plane_info(m_channel_map_name, detelement, plane);
    m_first_channel = static_cast<channel_t>(plane_info.min_channel);
    channel_t n_channels_on_plane = static_cast<channel_t>(plane_info.n_channels);

    // If we are in PD-VD and on the collection plane use 'effective' channel mapping for CRPs
    if (plane != 2) m_pdvd_map = false; // Effective mapping in PD-VD only for collection plane
    if (m_pdvd_map) {
      m_pdvd_eff_channel_mapper = std::make_unique<PDVDEffectiveChannelMap>(plane_info.min_channel, plane_info.n_channels);
      m_first_channel = m_pdvd_eff_channel_mapper->remapCollectionPlaneChannel(m_first_channel);
      m_last_channel = m_first_channel + m_pdvd_eff_channel_mapper->getNEffectiveChannels();
      m_chan_bin_length = m_pdvd_eff_channel_mapper->getNEffectiveChannels() / m_num_chanbins;
    } else {
      m_last_channel = m_first_channel + n_channels_on_plane;
      m_chan_bin_length = n_channels_on_plane / m_num_chanbins;
    }

    //PlaneInfo plane_info = get_plane_info(0, detelement, plane);

    //m_first_channel = static_cast<channel_t>(plane_info.min_channel);
    //channel_t n_channels_on_plane = static_cast<channel_t>(plane_info.n_channels);
    
    TLOG_DEBUG(TLVL_DEBUG_ALL) << "[TAM:PDBSMW] 1st Chan = " << m_first_channel << ", last Chan = " << m_last_channel << std::endl
      << "Number of channel bins = " << m_num_chanbins << ", and channel bin length = " << m_chan_bin_length;
    
    //std::cout << "IN TAMAKER: det, pl: " << detelement << "," << plane << " = min, n :" << m_first_channel << ", " << n_channels_on_plane << "\n";
    //std::cout << "1st Chan = " << m_first_channel << ", last Chan = " << m_last_channel << std::endl
    //  << "Number of bins = " << m_num_chanbins << ", and bin length = " << m_chan_bin_length << std::endl;
    
    return;
  } 

  // If the difference between the current TP's start time and the start of the window
  // is less than the specified window size, add the TP to the window.
  if((input_tp.time_start - m_current_window.time_start) < m_window_length){
    TLOG_DEBUG(TLVL_DEBUG_HIGH) << "[TAM:PDBSMW] Window not yet complete, adding the input_tp to the window.";
    m_current_window.add(input_tp);
  }
  // If the addition of the current TP to the window would make it longer
  // than the specified window length, don't add it
  // Instead go through a series of filters and eventually a XGBoost model to determine whether to create a TA
  else if (
      (m_current_window.time_start - m_last_pred_time) > m_bin_length && // check enough time has passed since last window
      m_current_window.tp_list.size() > 20 && // need enough TPs in window to bother
      m_current_window.adc_integral > m_adc_threshold && // set a low minimum threshold for the ADC integral sum
      (m_current_window.mean_adc_peak() / m_current_window.mean_tot()) > m_ratio_threshold && // mean peak / tot cut
      compute_treelite_classification() // XGBoost classifier 
      )
  {
    TLOG_DEBUG(TLVL_DEBUG_LOW) << "[TAM:PDBSMW] ADC integral in window is greater than specified threshold.";
    output_ta.push_back(construct_ta());
    TLOG_DEBUG(TLVL_DEBUG_HIGH) << "[TAM:PDBSMW] Resetting window with input_tp.";
    m_current_window.reset(input_tp);
  }
  // If it is not, move the window along.
  else{
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
    if (config.contains("channel_map_name")) m_channel_map_name = config["channel_map_name"];
    if (config.contains("num_time_bins")) m_num_timebins = config["num_time_bins"];
    if (config.contains("adc_threshold")) m_adc_threshold = config["adc_threshold"];
    if (config.contains("ratio_threshold")) {
      m_ratio_threshold = config["ratio_threshold"];
      m_ratio_threshold *= 0.01;
    }
    if (config.contains("window_length")) {
      m_window_length = config["window_length"];
      m_bin_length = static_cast<timestamp_t>(m_window_length / m_num_timebins);
    }
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
 
  TLOG_DEBUG(TLVL_DEBUG_ALL) << "[TAM:PDBSMW] Bin length is " << m_bin_length << " for a window of " << m_num_timebins <<
    " bins. ADC threshold across window set to " << m_adc_threshold;
  std::cout << "Bin length is " << m_bin_length << " for a window of " << m_num_timebins << 
    " bins. ADC threshold across window set to " << m_adc_threshold << std::endl;

  channelMap = dunedaq::detchannelmaps::make_map(m_channel_map_name);

  // If we are in PD-VD, set boolean to true to enable effective channel mapping
  if (m_channel_map_name == "PD2VDTPCChannelMap" || m_channel_map_name == "PD2VDBottomTPCChannelMap" ||
      m_channel_map_name == "PD2VDTopTPCChannelMap") {
    m_pdvd_map = true;
  } else { // else we are in PD-HD and we use true channel mapping
    m_pdvd_map = false;
  }

  m_compiled_model_interface = std::make_unique<CompiledModelInterface>(nbatch);

  const size_t num_feature = m_compiled_model_interface->GetNumFeatures();

  flat_batched_inputs.resize(num_feature);

  m_num_chanbins = num_feature / m_num_timebins;

  flat_batched_Entries.clear();
  for (size_t i = 0; i < num_feature; ++i) {
    union Entry zero;
    zero.fvalue = 0.0;
    flat_batched_Entries.emplace_back(zero);
  }

}

TriggerActivityMakerProtoDUNEBSMWindow::~TriggerActivityMakerProtoDUNEBSMWindow() {
  // Treelite smart ptr should clean itself up
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
      flat_batched_inputs, m_bin_length,
      m_num_chanbins, m_first_channel,
      m_chan_bin_length, m_num_timebins,
      m_pdvd_eff_channel_mapper, m_pdvd_map
      );

  m_current_window.fill_entry_window(flat_batched_Entries, flat_batched_inputs);

  std::vector<float> result(nbatch, 0.0f);

  m_compiled_model_interface->Predict(flat_batched_Entries.data(), result.data());

  return m_compiled_model_interface->Classify(result.data(), m_bdt_threshold);
}

// Register algo in TA Factory
REGISTER_TRIGGER_ACTIVITY_MAKER(TRACE_NAME, TriggerActivityMakerProtoDUNEBSMWindow)
