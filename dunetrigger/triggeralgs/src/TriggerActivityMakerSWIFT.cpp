#include "dunetrigger/triggeralgs/include/triggeralgs/SWIFT/TriggerActivityMakerSWIFT.hpp"

#include "TRACE/trace.h"
#define TRACE_NAME "TriggerActivityMakerSWIFTPlugin"

#include <cassert>
#include <iostream>
#include <cmath> 

namespace triggeralgs {

  // configuration
  void TriggerActivityMakerSWIFT::configure(const nlohmann::json& config)
  {
    //window settings
    if (config.contains("window_length"))
      m_window_length = config["window_length"];

    if (config.contains("inspect_energy_threshold"))
      m_inspect_energy_threshold = config["inspect_energy_threshold"];

    if (config.contains("accept_energy_threshold"))
      m_accept_energy_threshold = config["accept_energy_threshold"];

    //tp filtering settings
    if (config.contains("min_adc_peak"))
      m_min_adc_peak = config["min_adc_peak"];

    if (config.contains("min_samples_over_threshold"))
      m_min_samples_over_threshold = config["min_samples_over_threshold"]; 
    
    //clustering 
    if (config.contains("cm_per_tick"))
      m_cm_per_tick = config["cm_per_tick"];

    if (config.contains("wire_pitch"))
      m_wire_pitch = config["wire_pitch"];

    if (config.contains("min_samples"))
      m_db_min_samples = config["min_samples"];

    if (config.contains("epsilon")) // NN search radius 
      m_db_eps = config["epsilon"];

    if (config.contains("cluster_energy_cut"))
      m_cluster_energy_cut = config["cluster_energy_cut"];

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

  /*
  void TriggerActivityMakerSWIFT::close_window(std::vector<TriggerActivity>& output_tas) {
    int flag = 0;

    if (m_window_energy > m_accept_energy_threshold)
      flag = 2;
    else if (m_window_energy > m_inspect_energy_threshold)
      flag = 1;

    // Emit TA only if window passes thresholds
    if (flag > 0 && !m_current_ta.inputs.empty()) {
      set_ta_attributes();
      output_tas.push_back(m_current_ta);
    }

  }
  */

  //HERE 

  void TriggerActivityMakerSWIFT::close_window(std::vector<TriggerActivity>& output_tas) {
    int flag = 0;

    if (m_window_energy > m_accept_energy_threshold)
      flag = 2;
    else if (m_window_energy > m_inspect_energy_threshold)
      flag = 1;
    
    if (flag == 1) { // window requires inspection
      uint64_t max_cluster_energy = extract_dominant_cluster_energy(m_current_ta.inputs, m_db_eps, m_db_min_samples);

      if (max_cluster_energy > m_cluster_energy_cut)
        flag = 2; // upgrade window
    }

    // Emit TA only if window passes thresholds
    if (flag ==2  && !m_current_ta.inputs.empty()) {
      set_ta_attributes();
      output_tas.push_back(m_current_ta);
    }

  }
  
  //clustering - keeping only necessary for info for now, but might want to consider storing other properties & clusters themselves, not just max eng.   
  uint64_t TriggerActivityMakerSWIFT::extract_dominant_cluster_energy(const std::vector<TriggerPrimitive>& tps, float eps, int min_samples){
    
    // translate TriggerPrimitives -> Points
     struct Point {
       float z;   // channel -> z [cm]
       float t;   // time_start -> effective drift time i.e. x [cm]
       uint64_t adc;
     };

     std::vector<Point> points;
     points.reserve(tps.size());
     for (const auto& tp : tps) {
       points.push_back({tp.channel * m_wire_pitch,
	     (tp.time_start - m_window_start) * m_cm_per_tick,
	     tp.adc_integral});
     }


    //cluster labelling
    const int N = points.size();
    int cluster_id = 0;

    std::vector<int> labels(N, -1); // Points with label -1 = noise
    std::vector<uint8_t> visited(N, 0);

    std::vector<float> cluster_sums;       // store cluster energies

    float eps2 = eps * eps; //search radius 

    //nearest neighbour searches + clustering 
    for (int i = 0; i < N; ++i) {
      if (visited[i]) continue;
      visited[i] = 1;

      // find neighbors of i
      std::vector<int> neigh;
      const float zi = points[i].z;
      const float ti = points[i].t;

      for (int j = 0; j < N; ++j) {
	float dz = points[j].z - zi;
	float dt = points[j].t - ti;
	if (dz*dz + dt*dt <= eps2) {
	  neigh.push_back(j);
	}
      }
      
      //valid cluster only if min. number of points is reached
      if (neigh.size() < static_cast<size_t>(min_samples)) {
	labels[i] = -1;  // noise
	continue;
      }

      // starting new cluster
      labels[i] = cluster_id;
      uint64_t cluster_energy = points[i].adc;

      size_t k = 0;
      while (k < neigh.size()) {
	int j = neigh[k];

	if (!visited[j]) {
	  visited[j] = 1;

	  // expand neighbors of j
	  std::vector<int> neigh2;
	  const float zj = points[j].z;
	  const float tj = points[j].t;

	  for (int m = 0; m < N; ++m) {
	    float dz = points[m].z - zj;
	    float dt = points[m].t - tj;
	    if (dz*dz + dt*dt <= eps2) {
	      neigh2.push_back(m);
	    }
	  }

	  if (neigh2.size() >= static_cast<size_t>(min_samples)) {
	    for (int m : neigh2) {
	      if (std::find(neigh.begin(), neigh.end(), m) == neigh.end()) {
		neigh.push_back(m);
	      }
	    }
	  }
	}

	if (labels[j] == -1) {
	  labels[j] = cluster_id;
	  cluster_energy += points[j].adc;
	} else if (labels[j] == cluster_id) {
	  cluster_energy += points[j].adc;
	}

	++k;
      }

      cluster_sums.push_back(cluster_energy);
      ++cluster_id;
    }

    if (cluster_sums.empty()) return 0;

    return *std::max_element(cluster_sums.begin(), cluster_sums.end());
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
