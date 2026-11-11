#include "dunetrigger/experimental/tpctriggeralgs/include/tpctriggeralgs/SWIFT/TriggerActivityROPMakerSWIFT.hpp"

#include "dunetrigger/experimental/tpctriggeralgs/include/tpctriggeralgs/ConfigUtils.hpp"

#include "TRACE/trace.h"
#define TRACE_NAME "TriggerActivityROPMakerSWIFTPlugin"

#include <iostream>
#include <cmath>
#include <stdexcept>

namespace tpctriggeralgs {

  // configuration
  void TriggerActivityROPMakerSWIFT::configure(const nlohmann::json& config)
  {
    check_known_keys(config,
                     { "window_length", "inspect_energy_threshold_sadc", "accept_energy_threshold_sadc",
                       "min_adc_peak", "min_samples_over_threshold",
                       "cm_per_tick", "wire_pitch", "min_samples", "epsilon", "cluster_energy_cut_sadc" },
                     TRACE_NAME);

    // All parameters are required: no defaults
    //window settings
    m_window_length                 = get_required<uint64_t>(config, "window_length", TRACE_NAME); // in DTS ticks
    m_inspect_energy_threshold_sadc = get_required<uint64_t>(config, "inspect_energy_threshold_sadc", TRACE_NAME);
    m_accept_energy_threshold_sadc  = get_required<uint64_t>(config, "accept_energy_threshold_sadc", TRACE_NAME);

    //tp filtering settings
    m_min_adc_peak               = get_required<uint16_t>(config, "min_adc_peak", TRACE_NAME); //ADC
    m_min_samples_over_threshold = get_required<uint16_t>(config, "min_samples_over_threshold", TRACE_NAME); // in DTS ticks while still using TPv1 FIXME

    //clustering
    //these depend on detector properties
    m_cm_per_tick              = get_required<float>(config, "cm_per_tick", TRACE_NAME); // DTS tick period [us/tick] * drift velocity [cm/us]
    m_wire_pitch               = get_required<float>(config, "wire_pitch", TRACE_NAME); // cm
    m_db_min_samples           = get_required<int>(config, "min_samples", TRACE_NAME); //min. number of TPs for valid cluster
    m_db_eps                   = get_required<float>(config, "epsilon", TRACE_NAME); //dbscan search radius in cm
    m_cluster_energy_cut_sadc  = get_required<uint64_t>(config, "cluster_energy_cut_sadc", TRACE_NAME); // min energy of dominant cluster eng. in window for acceptance

    // window_length is used as a divisor; must be checked at runtime (asserts are compiled out in prof builds)
    if (m_window_length == 0)
      throw std::invalid_argument(TRACE_NAME ": window_length must be > 0");
  }

  //TP refinement
  bool TriggerActivityROPMakerSWIFT::preprocess( const TriggerPrimitive& input_tp) const{
    //FIXME: change TOT -> SOT after updating to TP v2.
    return (input_tp.adc_peak > m_min_adc_peak) && (input_tp.time_over_threshold > m_min_samples_over_threshold);
  }

  // Reset window state
  void TriggerActivityROPMakerSWIFT::reset_window_state(uint64_t new_window_start) {
    m_window_start       = new_window_start;
    m_window_energy_sadc = 0;
    m_tp_count           = 0;
    m_current_ta.inputs.clear();  // reuse existing TA heap allocation
    m_current_ta.time_start = m_window_start;
  }

  //main function for binning TPs into fixed-size time windows
  // **assumes TPs are strictly time ordered**
  void TriggerActivityROPMakerSWIFT::operator()(const TriggerPrimitive& input_tp, std::vector<TriggerActivityROP>& output_tas)
  {
    // Apply TP filtering
    if (!preprocess(input_tp)) {
      return;
    }

    // If TP is valid, determine which window this TP belongs to
    const uint64_t tp_window_start = (input_tp.time_start / m_window_length) * m_window_length;

    // Initialise on first TP
    if (!m_initialised) {

      //create TA instance once at run start
       m_current_ta = TriggerActivityROP();

      //reset window state
      reset_window_state(tp_window_start);
      m_initialised = true;
    }

    //if TP belongs to past window (already closed), reject it (not ideal)
    if (tp_window_start < m_window_start){
      return;
    }

    //if TP belongs to future window, close existing TA and start a new one at current time
    if (tp_window_start > m_window_start){
      close_window(output_tas);
      reset_window_state(tp_window_start);
    }

    //If we got here, the TP belongs to the current window
    m_current_ta.inputs.push_back(input_tp);
    m_window_energy_sadc += input_tp.adc_integral;
    ++m_tp_count;
  }

  //function which gets called when the window is closed.
  //main window categorisation & TA generation happens here
  void TriggerActivityROPMakerSWIFT::close_window(std::vector<TriggerActivityROP>& output_tas) {

    if (m_current_ta.inputs.empty()) return;

    //Prompt window categorisatoin : immidiate accept, inspect, reject based on local energy in window
    WindowDecision decision;
    if (m_window_energy_sadc > m_accept_energy_threshold_sadc) decision = WindowDecision::kAccept;
    else if (m_window_energy_sadc > m_inspect_energy_threshold_sadc) decision = WindowDecision::kInspect;
    else  return; // Reject

    // cluster inspect cases
    if (decision == WindowDecision::kInspect) {
      const uint64_t max_cluster_energy_sadc =  extract_dominant_cluster_energy(m_current_ta.inputs, m_db_eps, m_db_min_samples);
      if (max_cluster_energy_sadc <= m_cluster_energy_cut_sadc) return; // reject window if it didn't pass inspection
    }

    // Emit TA: should only reach this step if dealing with immidiate or conditional accept cases.
    set_ta_attributes();
    output_tas.push_back(m_current_ta);
  }

  //Clustering function: density-based clustering in t-z, where both coordinates are expressed in cm.
  //returns only max eng. for now, as that's the param. based on which decision is made.
  //FIXME eventually want this step to return nclusrers, mean cluster eng. etc.
  uint64_t TriggerActivityROPMakerSWIFT::extract_dominant_cluster_energy(const std::vector<TriggerPrimitive>& tps, float eps, int min_samples){

    struct Point {
      float z;
      float t;
      uint64_t adc;
    };

    //since time and channel are in different units, express TPs as points in unified coordinate system (z,t)
    std::vector<Point> points;
    points.reserve(tps.size());
    for (const auto& tp : tps) {
      points.push_back({tp.channel * m_wire_pitch, (tp.time_start - m_window_start) * m_cm_per_tick, tp.adc_integral});
    }

    const size_t N = points.size(); //number of elements to cluster
    int cluster_id = 0;
    std::vector<int> labels(N, -1); // initialise all labels to noise for now
    std::vector<uint8_t> visited(N, 0);
    float eps2 = eps * eps; //clustering radius

    for (size_t i = 0; i < N; ++i) {
      if (visited[i]) continue;
      visited[i] = 1;

      std::vector<int> neigh;
      const float zi = points[i].z;
      const float ti = points[i].t;

      for (size_t j = 0; j < N; ++j) {
        float dz = points[j].z - zi;
        float dt = points[j].t - ti;
        if (dz*dz + dt*dt <= eps2) neigh.push_back(j);
      }

      if (neigh.size() < static_cast<size_t>(min_samples)) {
        labels[i] = -1;
        continue;
      }

      labels[i] = cluster_id;

      size_t k = 0;
      while (k < neigh.size()) {
        int j = neigh[k];

        if (!visited[j]) {
          visited[j] = 1;

          std::vector<int> neigh2;
          const float zj = points[j].z;
          const float tj = points[j].t;
          for (size_t m = 0; m < N; ++m) {
            float dz = points[m].z - zj;
            float dt = points[m].t - tj;
            if (dz*dz + dt*dt <= eps2) neigh2.push_back(m);
          }

          if (neigh2.size() >= static_cast<size_t>(min_samples)) {
            for (int m : neigh2) {
              if (std::find(neigh.begin(), neigh.end(), m) == neigh.end())
                neigh.push_back(m);
            }
          }
        }

        if (labels[j] == -1) labels[j] = cluster_id;
        ++k;
      }

      ++cluster_id;
    }

    // once clusters are formed, calculate the energies
    std::vector<uint64_t> cluster_sums(cluster_id, 0);
    for (size_t i = 0; i < N; ++i) {
      if (labels[i] >= 0) cluster_sums[labels[i]] += points[i].adc;
    }
    //return dominant cluster energy in window
    if (cluster_sums.empty()) return 0;
    return *std::max_element(cluster_sums.begin(), cluster_sums.end());
  }



  //TA feature extraction - FIXME most of these fields are somewhat useless
  void TriggerActivityROPMakerSWIFT::set_ta_attributes()
  {
    m_current_ta.time_start = m_window_start;
    m_current_ta.time_end   = m_window_start + m_window_length;
    m_current_ta.adc_integral = m_window_energy_sadc;

    const TriggerPrimitive& first_tp = m_current_ta.inputs.front();
    m_current_ta.detid = first_tp.detid;
    m_current_ta.type = static_cast<uint8_t>(dunedaq::trgdataformats::exp::TriggerActivityROPData::Type::kTPC);
    m_current_ta.algorithm = static_cast<uint8_t>(dunedaq::trgdataformats::exp::TriggerActivityROPData::Algorithm::kUnknown);


    tpctriggeralgs::channel_t min_ch = first_tp.channel;
    tpctriggeralgs::channel_t max_ch = first_tp.channel;

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


  void TriggerActivityROPMakerSWIFT::flush(timestamp_t /*until*/, std::vector<TriggerActivityROP>& output_tas)
  {
    if (!m_initialised) return;
    close_window(output_tas);
  }

  REGISTER_TRIGGER_ACTIVITY_ROP_MAKER(TRACE_NAME, TriggerActivityROPMakerSWIFT)

} // namespace tpctriggeralgs
