/**
 * @file TriggerActivityMakerAEAnomalyWindow.hpp
 *
 * This is part of the DUNE DAQ Application Framework, copyright 2021.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */

#ifndef TRIGGERALGS_AEWINDOW_TRIGGERACTIVITYMAKERAEWINDOW_HPP_
#define TRIGGERALGS_AEWINDOW_TRIGGERACTIVITYMAKERAEWINDOW_HPP_

#include "dunetrigger/triggeralgs/include/triggeralgs/TriggerActivityFactory.hpp"
#include "dunetrigger/triggeralgs/include/triggeralgs/Types.hpp"
#include "dunetrigger/triggeralgs/include/triggeralgs/AEAnomalyWindow/WindowBin.hpp"
#include "dunetrigger/triggeralgs/include/triggeralgs/AEAnomalyWindow/BinnedWindow.hpp"

#include "tensorflow/core/public/session.h"
#include "tensorflow/core/platform/env.h"

#include <fstream>
#include <vector>
#include <algorithm>

namespace triggeralgs {
class TriggerActivityMakerAEAnomalyWindow : public TriggerActivityMaker
{

public:
  void operator()(const TriggerPrimitive& input_tp, std::vector<TriggerActivity>& output_ta);
  
  void configure(const nlohmann::json &config);

private:

  // Information on the model
  struct ModelConfig {
    std::string modelName;
    float normMin;
    float normMax;
    std::string inputTensor;
    std::string outputTensor;
  };

  ModelConfig loadModelConfig(const std::string &filename);

  std::vector<float> minMaxNormalize(const std::vector<float>& data);
 
  float computeMSE(const std::vector<float>& input, const std::vector<float>& output);

  // function to provide output data through nn inference
  std::vector<float> run_inference(std::vector<float> &input);

  TriggerActivity construct_ta() const;

  BinnedWindow m_current_window;
  WindowBin m_current_bin;
  uint64_t m_primitive_count = 0;

  // Configurable parameters.
  uint32_t m_adc_threshold = 1200000;
  float m_mse_threshold = 0.0001;
  // now the length of the bin
  // 250 tick bins and 80 bins in window
  timestamp_t m_bin_length = 250;
  int nbins = 80;
  bool m_use_ae = true;
  // model and tensorflow session parameters for inference
  ModelConfig m_model_config;
  std::string m_model;
  tensorflow::GraphDef graph_def;
  tensorflow::Session *m_session;
  tensorflow::Status m_status;
  tensorflow::Tensor m_input_tensor;

  // Try adding some TH1Ds for checks
  //TH1D *hInput;
  //TH1D *hOutput;

};
} // namespace triggeralgs

#endif // TRIGGERALGS_AEWINDOW_TRIGGERACTIVITYMAKERAEWINDOW_HPP_
