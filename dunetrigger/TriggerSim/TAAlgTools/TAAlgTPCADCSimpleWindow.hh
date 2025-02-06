#ifndef DUNETRIGGER_TRIGGERSIM_TAALGTPCADCSIMPLEWINDOW_hh
#define DUNETRIGGER_TRIGGERSIM_TAALGTPCADCSIMPLEWINDOW_hh

#include "fhiclcpp/ParameterSet.h" 

#include "detdataformats/DetID.hpp"
#include "dunetrigger/TriggerSim/TAAlgTools/TAAlgTPCTool.hh"
#include "dunetrigger/TriggerSim/TAAlgTools/TPWindow.hh"

#include "dunetrigger/TriggerSim/Verbosity.hh"

namespace duneana {

  class TAAlgTPCADCSimpleWindow : public TAAlgTPCTool {

  public:
    explicit TAAlgTPCADCSimpleWindow(fhicl::ParameterSet const& ps)
      : adc_threshold_(ps.get<uint32_t>("adc_threshold", 100000))
      , window_length_(ps.get<dunedaq::trgdataformats::timestamp_t>("window_length", 10000))
      , verbosity_(ps.get<int>("verbosity", 0)) //use it for couts
    {}
    
    dunedaq::trgdataformats::TriggerActivityData construct_ta() const
    {
      
      dunedaq::trgdataformats::TriggerActivityData ta;
      
      dunedaq::trgdataformats::TriggerPrimitive last_tp = current_window_.inputs.back();

      ta.time_start = last_tp.time_start;
      ta.time_end = last_tp.time_start;
      ta.time_peak = last_tp.time_peak;
      ta.time_activity = last_tp.time_peak;
      ta.channel_start = last_tp.channel;
      ta.channel_end = last_tp.channel;
      ta.channel_peak = last_tp.channel;
      ta.adc_integral = current_window_.adc_integral;
      ta.adc_peak = last_tp.adc_peak;
      ta.detid = last_tp.detid;
      ta.type = dunedaq::trgdataformats::TriggerActivityData::Type::kTPC;
      ta.algorithm = dunedaq::trgdataformats::TriggerActivityData::Algorithm::kADCSimpleWindow;

      for (const auto& tp : current_window_.inputs) {
	ta.time_start = std::min(ta.time_start, tp.time_start);
	ta.time_end = std::max(ta.time_end, tp.time_start);
	ta.channel_start = std::min(ta.channel_start, tp.channel);
	ta.channel_end = std::max(ta.channel_end, tp.channel);
	if (tp.adc_peak > ta.adc_peak) {
	  ta.time_peak = tp.time_peak;
	  ta.adc_peak = tp.adc_peak;
	  ta.channel_peak = tp.channel;
	}
      }
      
      return ta;
    }

    void initialize() override
    {
      ta_current_ = TriggerActivity();
      current_window_ = TPWindow();
    }

    //process TPs one-by-one
    //if adc integral is greater than the specified threshold for a given TP window, construct a TA, and add it to the output TAs vector
    void process_tp(art::Ptr<dunedaq::trgdataformats::TriggerPrimitive> tp,
		    std::vector<TriggerActivity> & tas_out) override
    {

      dunedaq::trgdataformats::TriggerPrimitive input_tp = *tp;

      // For the first TP, reset the window object.
      if (current_window_.is_empty()) {
	current_window_.reset(input_tp);
	ta_current_.second.push_back(tp);

	if (verbosity_ >= Verbosity::kInfo) std::cout << " TP start Time: " << input_tp.time_start << ", TP ADC Sum: " << input_tp.adc_integral
				      << ", TP TOT: " << input_tp.time_over_threshold << ", TP ADC Peak: " << input_tp.adc_peak
				      << ", TP Offline Channel ID: " << input_tp.channel << "\n";
	return;
      } 
      
      // If the difference between the current TP's start time and the window's start time
      // is less than the specified window length, add the TP to the window
      if ((input_tp.time_start - current_window_.time_start) < window_length_) {
	current_window_.add(input_tp);
	ta_current_.second.push_back(tp);
      }
      // If the addition of the current TP to the window makes it longer
      // than the specified window length, don't add it but check whether the adc integral in
      // the existing window is above the specified threshold. If true, make a TA and start 
      // a fresh window with the current TP.
      else if (current_window_.adc_integral > adc_threshold_) {
	ta_current_.first = construct_ta();
	
	tas_out.push_back(ta_current_);
	initialize();

	ta_current_.second.push_back(tp);	
	current_window_.reset(input_tp);
      }
      // If false, move the window along
      else {
	current_window_.move(input_tp, window_length_);
	ta_current_.second.push_back(tp);
	
	// removing tps from ta_current_.second which are removed from current_window_ using move()
	art::PtrVector<dunedaq::trgdataformats::TriggerPrimitive> tmp_tp_vec = ta_current_.second;
	ta_current_.second.clear();
	for (const auto tp_window : current_window_.inputs) {
	  for (const auto& tp_tmp : tmp_tp_vec) {
	    if (tp_window.channel == (*tp_tmp).channel &&
		tp_window.time_start == (*tp_tmp).time_start && 
		tp_window.time_over_threshold == (*tp_tmp).time_over_threshold &&
		tp_window.time_peak == (*tp_tmp).time_peak &&
		tp_window.adc_integral == (*tp_tmp).adc_integral && 
		tp_window.adc_peak == (*tp_tmp).adc_peak) {
	      ta_current_.second.push_back(tp_tmp);
	      break;
	    }
	  }
	}
      }
      return;      
    }
    
    
  private:
    uint32_t adc_threshold_;
    dunedaq::trgdataformats::timestamp_t window_length_;
    const int verbosity_;
    
    TPWindow current_window_;
    TriggerActivity ta_current_;
  };
}

#endif
