/**
 * @file TriggerActivityMakerSWIFT.hpp
 *
 * This is part of the DUNE DAQ Application Framework, copyright 2021.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */

#ifndef TRIGGERALGS_SWIFT_TRIGGERACTIVITYMAKERSWIFT_HPP_
#define TRIGGERALGS_SWIFT_TRIGGERACTIVITYMAKERSWIFT_HPP_

#include "dunetrigger/triggeralgs/include/triggeralgs/TriggerActivityFactory.hpp"
#include <algorithm>

namespace triggeralgs {
  
  class TriggerActivityMakerSWIFT : public TriggerActivityMaker {
  public:
    void operator()(const TriggerPrimitive& input_tp, std::vector<TriggerActivity>& output_tas);
    void configure(const nlohmann::json& config);
    void set_ta_attributes();
    void flush(std::vector<TriggerActivity>& output_tas);// override;
    bool preprocess( const TriggerPrimitive& input_tp);

  private:
    TriggerActivity m_current_ta;

    // configuration
    uint64_t m_window_length = 32000; // fixed window size in DTS ticks (32 * 1k readout ticks @ 500ns/tick)
    uint64_t m_inspect_energy_threshold = 15000;
    uint64_t m_accept_energy_threshold = 55000;

    //TP pre-processing
    uint16_t m_min_adc_peak = 80;
    uint16_t m_min_samples_over_threshold = 8;

       
    //clustering 
    // FIX ME : add epsilon and min. samples

    //window state
    bool     m_initialised = false;
    uint64_t m_window_start = 0;
    uint64_t m_window_energy = 0; 
    uint16_t m_tp_count = 0; 

    const timestamp_t kMaxTime = 6000 * 32;  //  maximum simulation time otherwise the windows keep rolling forever 

    void close_window(std::vector<TriggerActivity>& output_tas); 

};

} // namespace triggeralgs

#endif // TRIGGERALGS_SWIFT_TRIGGERACTIVITYMAKERSWIFT_HPP_
