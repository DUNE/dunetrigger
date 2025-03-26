/**
 * @file TriggerActivityMakerAEAnomalyWindow.cpp
 *
 * This is part of the DUNE DAQ Application Framework, copyright 2021.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */

#include "dunetrigger/triggeralgs/include/triggeralgs/AEAnomalyWindow/TriggerActivityMakerAEAnomalyWindow.hpp"

#include "TRACE/trace.h"
#define TRACE_NAME "TriggerActivityMakerAEAnomalyWindowPlugin"

#include <vector>


using namespace triggeralgs;
using Logging::TLVL_DEBUG_ALL;
using Logging::TLVL_DEBUG_HIGH;
using Logging::TLVL_DEBUG_LOW;
using Logging::TLVL_IMPORTANT;

void
TriggerActivityMakerAEAnomalyWindow::operator()(const TriggerPrimitive& input_tp, std::vector<TriggerActivity>& output_ta)
{
  
  // The first time operator is called, reset
  // window object.
  if(m_current_bin.is_empty()){
    m_current_bin.reset(input_tp);
    m_current_window.resetwindow(m_current_bin);
    m_primitive_count++;
    return;
  }
  //std::cout << "Running TriggerActivityMakerAEAnomalyWindowPlugin algorithm." << std::endl;

  // If the difference between the current TP's start time and the start of the bin
  // is less than the specified bin width, add the TP to the bin.
  if((input_tp.time_start - m_current_bin.time_start) < m_bin_length){
    TLOG_DEBUG(TLVL_DEBUG_HIGH) << "[TAM:ADCSW] Bin not yet complete, adding the input_tp to the bin.";
    //std::cout << "[TAM:ADCSW] Bin not yet complete, adding the input_tp to the bin." << std::endl;
    m_current_bin.add(input_tp);
  } 
  // If the window bins have not been all filled yet, add another bin
  else if (m_current_window.bincount() < nbins) {
    TLOG_DEBUG(TLVL_DEBUG_HIGH) << "[TAM:ADCSW] Window not yet complete, adding bin to the window.";
    std::cout << "[TAM:ADCSW] Window not yet complete at size " << m_current_window.bincount() << " out of " << nbins << ", adding bin to the window." << std::endl;
    m_current_window.addbin(m_current_bin);
    // Bin added - remember to reset the bin to start again
    m_current_bin.reset(input_tp);
  }
  // If the addition of the current TP to the window would make it longer
  // than the specified window length, don't add it but check whether the sum of all adc in
  // the existing window is above the specified threshold. If it is, make a TA and start 
  // a fresh window with the current TP.
  // This logic will be replaced with AE calculation
/*  else if(m_current_window.sumadc() > m_adc_threshold){
    TLOG_DEBUG(TLVL_DEBUG_LOW) << "[TAM:ADCSW] ADC integral in window is greater than specified threshold.";
    std::cout << "[TAM:ADCSW] ADC integral in window is greater than specified threshold." << std::endl;
    output_ta.push_back(construct_ta());
    TLOG_DEBUG(TLVL_DEBUG_HIGH) << "[TAM:ADCSW] Resetting window with input_tp.";
    std::cout << "[TAM:ADCSW] Resetting window with input_tp." << std::endl;
    m_current_bin.reset(input_tp);
    m_current_window.resetwindow(m_current_bin);
    std::cout << "Constructed TA with AE Window alg!" << std::endl;
  }
*/
  else if(computeMSE(m_current_window.ae_input, run_inference(m_current_window.ae_input)) > m_mse_threshold){
    std::cout << "mse calculated = " << computeMSE(m_current_window.ae_input, run_inference(m_current_window.ae_input)) 
      << ", and threshold = " << m_mse_threshold << std::endl;
    TLOG_DEBUG(TLVL_DEBUG_LOW) << "[TAM:ADCSW] MSE from AE reconstruction in window is greater than specified threshold.";
    std::cout << "[TAM:ADCSW] MSE from AE reconstruction in window is greater than specified threshold." << std::endl;
    output_ta.push_back(construct_ta());
    TLOG_DEBUG(TLVL_DEBUG_HIGH) << "[TAM:ADCSW] Resetting window with input_tp.";
    std::cout << "[TAM:ADCSW] Resetting window with input_tp." << std::endl;
    m_current_bin.reset(input_tp);
    m_current_window.resetwindow(m_current_bin);
    std::cout << "Constructed TA with AE Window alg!" << std::endl;
  }
  
  // If it is not, move the window along by removing the front bin and adding a new one to the back
  else{
    TLOG_DEBUG(TLVL_DEBUG_ALL) << "[TAM:ADCSW] Window is at required length but threshold not met, shifting window along by 1 window bin.";
    std::cout << "[TAM:ADCSW] Window is at required length but threshold not met, shifting window along by 1 window bin." << std::endl;
    //std::cout << "mse calculated = " << m_current_window.computeMSE(m_current_window.ae_input, run_inference(m_current_window.ae_input)) 
    //  << ", and threshold = " << m_mse_threshold << std::endl;
    //m_current_window.move(input_tp, m_bin_length);
    m_current_window.movebin(m_current_bin);
    // Do I then need to reset the current bin?
    // Probably... I have used the current bin and added it, now I need to start a new one
    m_current_bin.reset(input_tp);
  }
  
  //TLOG_DEBUG(TLVL_DEBUG_ALL) << "[TAM:ADCSW] " << m_current_window;

  m_primitive_count++;

  return;
}

void
TriggerActivityMakerAEAnomalyWindow::configure(const nlohmann::json &config)
{
  //FIXME use some schema here
  if (config.is_object()){
    if (config.contains("bin_length")) m_bin_length = config["bin_length"];
    if (config.contains("adc_threshold")) m_adc_threshold = config["adc_threshold"];
    if (config.contains("mse_threshold")) {
      uint64_t int_mse_threshold = config["mse_threshold"];
      m_mse_threshold = static_cast<float>(int_mse_threshold * 1E-3);
    }
    //if (config.contains("use_ae")) m_use_ae = config["use_ae"];
  }
  else{
    TLOG_DEBUG(TLVL_IMPORTANT) << "[TAM:ADCSW] The DEFAULT values of window_length and adc_threshold are being used.";
  }
  TLOG_DEBUG(TLVL_IMPORTANT) << "[TAM:ADCSW] If the total ADC of trigger primitives with times within a "
                         << m_bin_length << " tick time window is above " << m_adc_threshold << " counts, a trigger will be issued.";
  std::cout << "bin length is " << m_bin_length << " for a window if " << nbins << " bins. ADC threshold across window set to " << m_adc_threshold << std::endl;
  
  m_model = "/exp/dune/app/users/chasnip/CERN_Fellowship/dunetrigger_aeta/ae_window_model_1dcnn_dnnlayer_graph.pb";
  std::cout << "model loading: " << m_model << std::endl;
  
  // Create new session
  const tensorflow::SessionOptions options = tensorflow::SessionOptions();
  m_status = tensorflow::NewSession(options, &m_session);
  if (!m_status.ok()) {
    std::cerr << "Error creating session: " << m_status.ToString() << std::endl;
    return;
  }

  // Load model
  m_status = tensorflow::ReadBinaryProto(
          tensorflow::Env::Default(), m_model, &graph_def);
  if (!m_status.ok()) {
    std::cerr << "Error loading model: " << m_status.ToString() << std::endl;
    return;
  }
  
  // Add graph to session
  m_status = m_session->Create(graph_def);
  if (!m_status.ok()) {
    std::cerr << "Error adding graph: " << m_status.ToString() << std::endl;
    return;
  }

  // Create an input tensor with shape {1, 80, 1}
  m_input_tensor = tensorflow::Tensor(tensorflow::DT_FLOAT, tensorflow::TensorShape({1, nbins, 1}));
}

TriggerActivity
TriggerActivityMakerAEAnomalyWindow::construct_ta() const
{
  TLOG_DEBUG(TLVL_DEBUG_LOW) << "[TAM:ADCSW] I am constructing a trigger activity!";
  //TLOG_DEBUG(TRACE_NAME) << m_current_window;

  TriggerPrimitive latest_tp_in_window = m_current_bin.tp_list.back();
  // The time_peak, time_activity, channel_* and adc_peak fields of this TA are irrelevent
  // for the purpose of this trigger alg.
  TriggerActivity ta;
  ta.time_start = m_current_window.window_time_start;
  ta.time_end = latest_tp_in_window.time_start + latest_tp_in_window.time_over_threshold;
  ta.time_peak = latest_tp_in_window.time_peak;
  ta.time_activity = latest_tp_in_window.time_peak;
  ta.channel_start = latest_tp_in_window.channel;
  ta.channel_end = latest_tp_in_window.channel;
  ta.channel_peak = latest_tp_in_window.channel;
  ta.adc_integral = m_current_window.sumadc();
  ta.adc_peak = latest_tp_in_window.adc_peak;
  ta.detid = latest_tp_in_window.detid;
  ta.type = TriggerActivity::Type::kTPC;
  //ta.algorithm = TriggerActivity::Algorithm::kAEAnomalyWindow;
  ta.algorithm = TriggerActivity::Algorithm::kUnknown;
  // TODO: save the all the TPs, not just the ones in the latest bin
  //ta.inputs = m_current_window.tp_list;
  ta.inputs = m_current_bin.tp_list;
  return ta;
}

// Register algo in TA Factory
REGISTER_TRIGGER_ACTIVITY_MAKER(TRACE_NAME, TriggerActivityMakerAEAnomalyWindow)
