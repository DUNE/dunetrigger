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
        ae_input.push_back(input_bin.adc_integral);
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
      };

      float sumadc() const {
        float sum = std::accumulate(ae_input.begin(), ae_input.end(), 0);
        return sum;
      };

      int bincount() const {
        return (int)ae_input.size();
      };

      timestamp_t window_time_start;
      std::vector<WindowBin> tp_window_bins;
      std::vector<float> ae_input;
  };

  TriggerActivity construct_ta() const;

  BinnedWindow m_current_window;
  WindowBin m_current_bin;
  uint64_t m_primitive_count = 0;

  // Configurable parameters.
  uint32_t m_adc_threshold = 1200000;
  // now the length of the bin
  // 250 tick bins and 80 bins in window
  timestamp_t m_bin_length = 250;
  int nbins = 80;
};
} // namespace triggeralgs

#endif // TRIGGERALGS_ADCSIMPLEWINDOW_TRIGGERACTIVITYMAKERADCSIMPLEWINDOW_HPP_
