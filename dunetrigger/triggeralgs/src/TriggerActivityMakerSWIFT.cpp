#include "dunetrigger/triggeralgs/include/triggeralgs/SWIFT/TriggerActivityMakerSWIFT.hpp"

#include "TRACE/trace.h"
#define TRACE_NAME "TriggerActivityMakerSWIFTPlugin"

#include <cassert>
#include <iostream>

namespace triggeralgs {

  // configuration
  void TriggerActivityMakerSWIFT::configure(const nlohmann::json& config)
  {
    if (config.contains("window_length"))
      m_window_length = config["window_length"];

    if (config.contains("inspect_energy_threshold"))
      m_inspect_energy_threshold = config["inspect_energy_threshold"];

    if (config.contains("accept_energy_threshold"))
      m_accept_energy_threshold = config["accept_energy_threshold"];

    if (config.contains("min_adc_peak"))
      m_min_adc_peak = config["min_adc_peak"];

    if (config.contains("min_samples_over_threshold"))
      m_min_samples_over_threshold = config["min_samples_over_threshold"]; 
    
    assert(m_window_length > 0);
  }

  //TP refinement 
  bool TriggerActivityMakerSWIFT::preprocess( const TriggerPrimitive& input_tp){
    if ((input_tp.adc_peak < m_min_adc_peak) & (input_tp.time_over_threshold < m_min_samples_over_threshold * 32 )){ //FIXME - old TP format still using time over threshold, should be samples 
      return false;
    }
    return true;
  }


  //main window categorisation 
  void TriggerActivityMakerSWIFT::operator()(const TriggerPrimitive& input_tp, std::vector<TriggerActivity>& output_tas)
  {
    // Initialise global window clock
    if (!m_initialised) {
      // absolute window alignment across APA planes (wouldn't be the case if we used first tp time). This should probably be set to run start time FIXME
      m_window_start = 0;  
      m_initialised = true;

      //reset window state
      m_window_energy = 0;
      m_tp_count = 0;
      m_current_ta = TriggerActivity();
      m_current_ta.time_start = m_window_start;

    }

    // Apply TP filtering
    if (!preprocess(input_tp)) {
      return;
    }

    // Determine which window this TP belongs to
    const uint64_t tp_window_start = (input_tp.time_start / m_window_length) * m_window_length;
    
    // Advance windows until TP fits
    while (m_window_start < tp_window_start) {
      
      // Close current window (even if empty)
      close_window(output_tas);
      
      // Reset state for next window
      m_window_energy = 0;
      m_tp_count = 0;
      m_current_ta = TriggerActivity();
      m_current_ta.time_start = m_window_start;

      // Advance "global" clock
      m_window_start += m_window_length;

      // Safety guard so we don't run forever in simulations
      if (m_window_start >= kMaxTime) {
	return;
      }
    }

    // If TP belongs to the current window, add it. 
    m_current_ta.inputs.push_back(input_tp);
    m_window_energy += input_tp.adc_integral;
    ++m_tp_count;
  }


  void TriggerActivityMakerSWIFT::close_window(std::vector<TriggerActivity>& output_tas) {
    int flag = 0;

    if (m_window_energy > m_accept_energy_threshold)
      flag = 2;
    else if (m_window_energy > m_inspect_energy_threshold)
      flag = 1;

    // Emit TA only if window passes thresholds
    if (flag > 0 && !m_current_ta.inputs.empty()) {
      
      set_ta_attributes(); // save info about the TA
      //std::cout << flag << ", " << m_window_energy << ", " << m_tp_count << ", " << m_window_start << std::endl; 
      output_tas.push_back(m_current_ta);
    }

  }

  //TA feature extraction - FIXME most of these fields are somewhat useless 
  void TriggerActivityMakerSWIFT::set_ta_attributes()
  {
    m_current_ta.time_start = m_window_start;
    m_current_ta.time_end   = m_window_start + m_window_length;
    m_current_ta.adc_integral = m_window_energy;

    const TriggerPrimitive& first_tp = m_current_ta.inputs.front();
    m_current_ta.detid = first_tp.detid;
    m_current_ta.type = TriggerActivity::Type::kTPC;
    m_current_ta.algorithm = TriggerActivity::Algorithm::kUnknown;

   
    dunedaq::trgdataformats::channel_t min_ch = first_tp.channel;
    dunedaq::trgdataformats::channel_t max_ch = first_tp.channel;

    // Peak quantities
    m_current_ta.adc_peak = 0;
    for (const auto& tp : m_current_ta.inputs) {
      if (tp.channel < min_ch) min_ch = tp.channel;
      if (tp.channel > max_ch) max_ch = tp.channel;

      if (tp.adc_peak > m_current_ta.adc_peak) {
	m_current_ta.adc_peak = tp.adc_peak;
	m_current_ta.channel_peak = tp.channel;
	m_current_ta.time_peak = tp.time_peak;
      }
    }

    m_current_ta.channel_start = min_ch;
    m_current_ta.channel_end   = max_ch;
    m_current_ta.time_activity = m_current_ta.time_peak;

  }



  void TriggerActivityMakerSWIFT::flush(std::vector<TriggerActivity>& output_tas)
  {
    if (!m_initialised)
      return;

    close_window(output_tas);
  }

  REGISTER_TRIGGER_ACTIVITY_MAKER(TRACE_NAME, TriggerActivityMakerSWIFT)

} // namespace triggeralgs
