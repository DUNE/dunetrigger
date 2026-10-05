#ifndef DUNETRIGGER_TRIGGERSIM_TPALGTPCEXAMPLE_hh
#define DUNETRIGGER_TRIGGERSIM_TPALGTPCEXAMPLE_hh

#include "fhiclcpp/ParameterSet.h" 

#include "dunetrigger/TriggerSim/TPAlgTools/TPAlgTPCTool.hh"

namespace dunetrigger {

  class TPAlgTPCExample : public TPAlgTPCTool {

  public:
    explicit TPAlgTPCExample(fhicl::ParameterSet const& ps) :
      verbosity_(ps.get<int>("verbosity",0))
    {}

    void process_waveform(std::vector<short> const& adcs, 
                          triggeralgs::channel_t const channel,
                          triggeralgs::detid_t const detid,
                          triggeralgs::timestamp_t const start_time,
                          std::vector<triggeralgs::TriggerPrimitive> & tps_out) 
    {
      
      //do processing of waveform here to form a TP and add it to the output collection
      
      triggeralgs::TriggerPrimitive this_tp;
      
      this_tp.channel = channel;
      this_tp.detid = detid;
      
      this_tp.type =      triggeralgs::TriggerPrimitive::Type::kTPC;
      this_tp.algorithm = triggeralgs::TriggerPrimitive::Algorithm::kUnknown;
      
      this_tp.flag = 0;
      
      this_tp.time_start = start_time;
      this_tp.time_over_threshold = start_time + adcs.size()*this->ADC_SAMPLING_RATE_IN_DTS;
      this_tp.time_peak = adcs.size()*this->ADC_SAMPLING_RATE_IN_DTS / 2;
      
      this_tp.adc_integral = 0;
      this_tp.adc_peak = 0;
      
      tps_out.push_back(this_tp);
      
    }
    
  private:
    int verbosity_;
    
  };
  
}

#endif
