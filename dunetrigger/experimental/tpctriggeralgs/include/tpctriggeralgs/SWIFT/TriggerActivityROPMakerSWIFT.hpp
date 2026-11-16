/*
 * @file TriggerActivityROPMakerSWIFT.hpp
 *
 * This is part of the DUNE DAQ Application Framework, copyright 2021.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */

#ifndef tpctriggeralgs_SWIFT_TRIGGERACTIVITYROPMAKERSWIFT_HPP_
#define tpctriggeralgs_SWIFT_TRIGGERACTIVITYROPMAKERSWIFT_HPP_

#include "dunetrigger/experimental/tpctriggeralgs/include/tpctriggeralgs/TriggerActivityROPFactory.hpp"
#include <algorithm>

namespace tpctriggeralgs {

  class TriggerActivityROPMakerSWIFT : public TriggerActivityROPMaker {
  public:
    void operator()(const TriggerPrimitive& input_tp, std::vector<TriggerActivityROP>& output_tas) override;
    void configure(const nlohmann::json& config) override;
    void flush(timestamp_t until, std::vector<TriggerActivityROP>& output_tas) override;
    bool preprocess( const TriggerPrimitive& input_tp) const;

  private:
    TriggerActivityROP m_current_ta;
    
    // -- Algortihm configuration --
    // static variables which get set once at run start & kept fixed afterwards

    //Window settings
    uint64_t m_window_offset = 0; // in DTS ticks, used to shift the window start time relative to the first TP in the run. This is useful for aligning the window with the beam spill.
    uint64_t m_window_length;
    uint64_t m_inspect_energy_threshold_sadc; 
    uint64_t m_accept_energy_threshold_sadc;

    //TP pre-processing
    uint16_t m_min_adc_peak;
    uint16_t m_min_samples_over_threshold;

    //Clustering settings
    float m_cm_per_tick;
    float m_wire_pitch;
    int   m_db_min_samples;
    float m_db_eps;
    uint64_t m_cluster_energy_cut_sadc;


    // Dynamic vaiables which get updated during running to keep track of alg state. These are non-configurable

    //window state
    bool     m_initialised = false;
    uint64_t m_window_start = 0;
    uint64_t m_window_id = 0;
    uint64_t m_window_energy_sadc = 0;
    uint16_t m_window_tp_count = 0;
    std::vector<TriggerPrimitive> m_window_tps;

    //Additional SWIFT-specific functions and classes
    enum class WindowDecision { kReject, kInspect, kAccept };
    void close_window(std::vector<TriggerActivityROP>& output_tas);
    void reset_window_state(uint64_t new_window_start, uint64_t new_window_id);
    // uint64_t extract_dominant_cluster_energy(const std::vector<TriggerPrimitive>& tps, float eps, int min_samples);

    struct ClusterLabels {
      size_t num;
      std::vector<uint16_t> labels;
    };

    struct DominantClusterSADC {
      uint64_t sadc;
      size_t label_id;
    };

    struct Momenta {
      float ch_mean;
      float ch_std;
      float t_mean;
      float t_std;
      float ch_wmean;
      float ch_wstd;
      float t_wmean;
      float t_wstd;
    };

    // Calculates the cluster labels for each TP using DBSCAN clustering in t-z space, where both coordinates are expressed in cm.
    ClusterLabels find_cluster_labels(const std::vector<TriggerPrimitive>& tps, float eps, int min_samples);
    // Calculates the energy of each cluster and returns the maximum cluster energy
    DominantClusterSADC calculate_dominant_cluster_sadc(const std::vector<TriggerPrimitive>& tps, const ClusterLabels& cluster_result);
    // Calculate the mean and standard deviation of the channel and time of the TPs in the cluster, as well as the weighted mean and standard deviation using the ADC integral as weights.
    Momenta calculate_momenta(const std::vector<TriggerPrimitive>& tps, const std::function<bool(size_t)>& selector);

    // 
    void set_ta_attributes(const Momenta& momenta);
  };
} // namespace tpctriggeralgs

#endif // tpctriggeralgs_SWIFT_TRIGGERACTIVITYROPMAKERSWIFT_HPP_
