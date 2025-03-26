/**
 * @file TriggerActivityMakerAEAnomalyWindow.hpp
 *
 * This is part of the DUNE DAQ Application Framework, copyright 2021.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */

#ifndef TRIGGERALGS_ADCSIMPLEWINDOW_TRIGGERACTIVITYMAKERADCSIMPLEWINDOW_HPP_
#define TRIGGERALGS_ADCSIMPLEWINDOW_TRIGGERACTIVITYMAKERADCSIMPLEWINDOW_HPP_

#include "dunetrigger/triggeralgs/include/triggeralgs/TriggerActivityFactory.hpp"
#include "dunetrigger/triggeralgs/include/triggeralgs/Types.hpp"

#include "tensorflow/core/public/session.h"
#include "tensorflow/core/platform/env.h"

#include <vector>
#include <algorithm>

namespace triggeralgs {
class TriggerActivityMakerAEAnomalyWindow : public TriggerActivityMaker
{

public:
  void operator()(const TriggerPrimitive& input_tp, std::vector<TriggerActivity>& output_ta);
  
  void configure(const nlohmann::json &config);

private:

  // Make a second class for a mini-window for each bin that goes into the AE input
  // Re-used from original Window class in ADCSimpleWindow
  // Represents a single bin of TPs in time that goes into a window that is an input to an AE
  class WindowBin {
    public:
      bool is_empty() const{
        return tp_list.empty();
      };
      void add(TriggerPrimitive const &input_tp){
        // Add the input TP's contribution to the total ADC and add it to
        // the TP list.
        adc_integral += input_tp.adc_integral;
        tp_list.push_back(input_tp);
      };
      void clear(){
        tp_list.clear();
      };
      void move(TriggerPrimitive const &input_tp, timestamp_t const &window_length){
        // Find all of the TPs in the window that need to be removed
        // if the input_tp is to be added and the size of the window
        // is to be conserved.
        // Substract those TPs' contribution from the total window ADC.
        uint32_t n_tps_to_erase = 0;
        for(auto tp : tp_list){
          if(!(input_tp.time_start-tp.time_start < window_length)){
            n_tps_to_erase++;
            adc_integral -= tp.adc_integral;
          }
          else break;
        }
        // Erase the TPs from the window.
        tp_list.erase(tp_list.begin(), tp_list.begin()+n_tps_to_erase);
        // Make the window start time the start time of what is now the
        // first TP.
        if(tp_list.size()!=0){
          time_start = tp_list.front().time_start;
          add(input_tp);
        }
        else reset(input_tp);
      };
      void reset(TriggerPrimitive const &input_tp){
        // Empty the TP list.
        tp_list.clear();
        // Set the start time of the window to be the start time of the 
        // input_tp.
        time_start = input_tp.time_start;
        // Start the total ADC integral.
        adc_integral = input_tp.adc_integral;
        // Add the input TP to the TP list.
        tp_list.push_back(input_tp);
      };
      friend std::ostream& operator<<(std::ostream& os, const WindowBin& window){
        if(window.is_empty()) os << "Window is empty!\n";
        else{
          os << "Window start: " << window.time_start << ", end: " << window.tp_list.back().time_start;
          os << ". Total of: " << window.adc_integral << " ADC counts with " << window.tp_list.size() << " TPs.\n"; 
        }
        return os;
      };

      timestamp_t time_start;
      uint32_t adc_integral;
      std::vector<TriggerPrimitive> tp_list;
  };

  // Full window of TPs that are binned as a function of time
  class BinnedWindow {
    public:
      void addbin(WindowBin const &input_bin) {
        tp_window_bins.push_back(input_bin);
        //ae_input.push_back(input_bin.adc_integral);
        // Want floats ready for normalisation
        ae_input.push_back(static_cast<float>(input_bin.adc_integral));
        // normalise ae_input
        //minMaxNormalize(ae_input);
      };
      void resetwindow(WindowBin const &input_bin) {
        window_time_start = input_bin.time_start;
        tp_window_bins.clear();
        ae_input.clear();
      };
      void movebin(WindowBin const &input_bin) {
        // Add the next bin
        // All bins of equal length in time - so just pop out the front bin of TPs
        tp_window_bins.erase(tp_window_bins.begin());
        tp_window_bins.push_back(input_bin);

        ae_input.erase(ae_input.begin());
        ae_input.push_back(input_bin.adc_integral);
        // normalise ae_input
        //minMaxNormalize(ae_input);
      };

      float sumadc() const {
        float sum = std::accumulate(ae_input.begin(), ae_input.end(), 0);
        return sum;
      };

      int bincount() const {
        return (int)ae_input.size();
      };

      /*
      float computeMSE(const std::vector<float>& input, const std::vector<float>& output) {
        if (input.size() != output.size()) {
          std::cerr << "Error: Input and output vectors must have the same size." << std::endl;
          return -1.0;
        }

        float mse = 0.0;
        for (size_t i = 0; i < input.size(); i++) {
          float diff = input[i] - output[i];
          mse += diff * diff;
        }
        return mse / input.size();
      };
      */
      /*
      void minMaxNormalize(std::vector<float>& data) {
        if (data.empty()) return;
        // Find the min and max in the vector
        auto [min_it, max_it] = std::minmax_element(data.begin(), data.end());
        float min_val = *min_it;
        float max_val = *max_it;
        std::cout << "min = " << min_val << ", max = " << max_val << std::endl;
        // Prevent division by zero if all values are equal
        if (max_val == min_val) {
          std::fill(data.begin(), data.end(), 0.0f);
          return;
        }

        // Normalize each element: (x - min) / (max - min)
        for (float &x : data) {
          x = (x - min_val) / (max_val - min_val);
          std::cout << x << ", ";
        }
        std::cout << std::endl;
      }*/

      timestamp_t window_time_start;
      std::vector<WindowBin> tp_window_bins;
      std::vector<float> ae_input;
  };

  //void minMaxNormalize(std::vector<float>& data) {
  std::vector<float> minMaxNormalize(const std::vector<float>& data) {
    if (data.empty()) return {};

    // define return normalised vector
    std::vector<float> data_norm(data.size());

    // Find the min and max in the vector
    auto [min_it, max_it] = std::minmax_element(data.begin(), data.end());
    float min_val = *min_it;
    float max_val = *max_it;
    //std::cout << "min = " << min_val << ", max = " << max_val << std::endl;
    // Prevent division by zero if all values are equal
    if (max_val == min_val) {
      std::fill(data_norm.begin(), data_norm.end(), 0.0f);
      return data_norm;
    }

    // Normalize each element: (x - min) / (max - min)
    //for (float &x : data) {
    for (int i = 0; i < (int)data.size(); i++) {
      data_norm[i] = (data[i] - min_val) / (max_val - min_val);
      //x = (x - min_val) / (max_val - min_val);
      //std::cout << x << ", ";
    }
    //std::cout << std::endl;
    return data_norm;
  }
      
  float computeMSE(const std::vector<float>& input, const std::vector<float>& output) {
    if (input.size() != output.size()) {
      std::cerr << "Error: Input and output vectors must have the same size." << std::endl;
      return -1.0;
    }

    auto input_norm = minMaxNormalize(input);

    float mse = 0.0;
    for (size_t i = 0; i < input_norm.size(); i++) {
      float diff = input_norm[i] - output[i];
      mse += diff * diff;
    }
    return mse / input_norm.size();
  };

  // function to provide output data through nn inference
  std::vector<float> run_inference(std::vector<float> &input) {
    // Normalise input data
    auto input_norm = minMaxNormalize(input);
    
    for (size_t i = 0; i < input_norm.size(); i++) {
      m_input_tensor.flat<float>()(i) = input_norm[i];
    }
    std::vector<std::pair<std::string, tensorflow::Tensor>> t_inputs = {
      {"keras_tensor_36", m_input_tensor} };
    std::vector<tensorflow::Tensor> outputs;
    m_status = m_session->Run(t_inputs, {"StatefulPartitionedCall/Identity:0"}, {}, &outputs);

    if (!m_status.ok()) {
      std::cerr << "Error during inference: " << m_status.ToString() << std::endl;
    }

    //std::cout << "Output tensor type: " << outputs[0].dtype() << std::endl;
    //std::cout << "Output tensor shape: " << outputs[0].shape().DebugString() << std::endl;

    if (outputs.empty()) {
      std::cerr << "Error: No outputs were returned from the session." << std::endl;
    }

    if (outputs[0].dtype() != tensorflow::DT_FLOAT) {
      std::cerr << "Error: Output tensor is not of type float." << std::endl;
      return {};
    }

    auto output_tensor_map = outputs[0].flat<float>();
    std::vector<float> output_data(output_tensor_map.data(), output_tensor_map.data() + output_tensor_map.size());
    //std::cout << "INPUT DATA: \n";
    //for (const auto &i : input_norm) std::cout << i << ", ";
    //std::cout << "OUTPUT DATA: \n";
    //for (const auto &o : output_data) std::cout << o << ", ";
    //std::cout << std::endl;
    return output_data;
  };

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
  bool m_use_ae = false;
  std::string m_model;
  // tensorflow session for inference
  tensorflow::GraphDef graph_def;
  tensorflow::Session *m_session;
  tensorflow::Status m_status;
  tensorflow::Tensor m_input_tensor;

};
} // namespace triggeralgs

#endif // TRIGGERALGS_ADCSIMPLEWINDOW_TRIGGERACTIVITYMAKERADCSIMPLEWINDOW_HPP_
