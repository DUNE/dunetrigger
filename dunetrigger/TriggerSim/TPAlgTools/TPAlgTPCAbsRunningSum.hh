#ifndef DUNETRIGGER_TRIGGERSIM_ABSRUNNINGSUM_hh
#define DUNETRIGGER_TRIGGERSIM_ABSRUNNINGSUM_hh

#include "fhiclcpp/ParameterSet.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/WireReadout.h"

#include "dunetrigger/TriggerSim/TPAlgTools/TPAlgTPCTool.hh"
#include "dunetrigger/TriggerSim/Verbosity.hh"
#include <inttypes.h>
#include <limits>
#include <map>

namespace dunetrigger {

class TPAlgTPCAbsRunningSum : public TPAlgTPCTool {

public:
  explicit TPAlgTPCAbsRunningSum(fhicl::ParameterSet const &ps)
      : verbosity_(ps.get<int>("verbosity", 0)),
        accum_limit_(ps.get<int>("accum_limit", 10)),
        scale_factor_(ps.get<int>( "scale_factor", 5)), // AbsRS-specific parameter ensuring the scale of hit parameters remains in the same ballpark as raw ADCs 
        r_value_(ps.get<int>("r_value", 9)), // AbsRS-specific "memory" parameter.
        threshold_tpg_plane0_(ps.get<int16_t>("threshold_tpg_plane0")),
        threshold_tpg_plane1_(ps.get<int16_t>("threshold_tpg_plane1")),
        threshold_tpg_plane2_(ps.get<int16_t>("threshold_tpg_plane2")) {}

  void initialize_channel_state(dunedaq::trgdataformats::channel_t const &channel, std::vector<short> const &adcs) {

    // grab the geometry service
    art::ServiceHandle<geo::Geometry> geom;
    geo::WireReadoutGeom const &wireReadout = art::ServiceHandle<geo::WireReadout>()->Get();

    auto plane = wireReadout.ROPtoWirePlanes(wireReadout.ChannelToROP(channel)).at(0).Plane;

    if (plane == 0)
      threshold_ = threshold_tpg_plane0_;
    else if (plane == 1)
      threshold_ = threshold_tpg_plane1_;
    else if (plane == 2)
      threshold_ = threshold_tpg_plane2_;

    if (verbosity_ >= Verbosity::kInfo) {
      std::cout << "Channel: " << channel << ", ROP: " << plane
                << ", Threshold: " << threshold_ << std::endl;
    }

    // find the mode of the whole vector
    std::map<short, size_t> counts_per_value;
    for (auto const &adc : adcs)
      counts_per_value[adc] += 1;
    size_t max_counts = 0;
    for (auto it = counts_per_value.begin(); it != counts_per_value.end();++it) {
      if (it->second > max_counts) {
        pedestal_ = it->first;
        max_counts = it->second;
      }
    }

    accum_ = 0;

    // initialise the running sum/second pedestal variables for this waveform
    running_sum_ = 0;
    abs_pedestal_ = 0;
    abs_accum_ = 0;

    prev_was_over_ = false;
    hit_charge_ = 0;
    hit_tover_ = 0;
    hit_peak_adc_ = 0;
    hit_peak_time_ = 0;
  }

  void frugal_accum_update(const int16_t sample, int16_t &pedestal,int16_t &accum) {
    if (sample > pedestal)
      ++accum;
    if (sample < pedestal)
      --accum;

    if (accum > accum_limit_) {
      ++pedestal;
      accum = 0;
    }
    if (accum < -1 * accum_limit_) {
      --pedestal;
      accum = 0;
    }
  }

  // Absolute-value running sum
  void abs_run_sum(const int16_t sample) {
    // R-value should be a float in range [0,1], but we're working with integers
    // online. Divide using online division first, then scale up.
    const int16_t tmp_sum = r_value_ * avx2_divide(running_sum_, 10) + avx2_divide(std::abs(sample), 10) * scale_factor_;

    // Saturated additions. AbsRS should be strictly positive, so a negative tmp_sum indicates overflow.
    if (tmp_sum < 0)
      running_sum_ = std::numeric_limits<int16_t>::max();  // Set the sum to the saturation limit.
    else
      running_sum_ = tmp_sum;
  }


  void process_waveform(std::vector<short> const &adcs, dunedaq::trgdataformats::channel_t const channel,
			dunedaq::trgdataformats::detid_t const detid,  dunedaq::trgdataformats::timestamp_t const start_time,
			std::vector<dunedaq::trgdataformats::TriggerPrimitive> &tps_out) {

    // setup a TP and initialize it with the common things for this
    // algorithm/channel
    dunedaq::trgdataformats::TriggerPrimitive this_tp;

    this_tp.channel = channel;
    this_tp.detid = detid;
    this_tp.type = dunedaq::trgdataformats::TriggerPrimitive::Type::kTPC;
    this_tp.algorithm = dunedaq::trgdataformats::TriggerPrimitive::Algorithm::kAbsRunningSum;
    this_tp.flag = 0;

    // for this channel, reinitialize the channel state variables
    initialize_channel_state(channel, adcs);


    // Dry-run Phase: Process the first ~300 samples (no hit detection, only pedestal stabilization)
    // once that's done, process the waveform normally from the top, otherwise we get a lot of noise TPs before the baseline plateaus
    // this wouldn't be necessary in the online system where we only get 'dead time' at the beginning of the run, and then we have continuous readout 
    // but it becomes noticable in the simulation where the waveform has to be initialised over and over for each event  
    for (size_t jj = 0; jj < 300 ; ++jj) {
      int16_t sample = adcs[jj];

      // Update the pedestal estimate
      frugal_accum_update(sample, pedestal_, accum_);
      sample -= pedestal_;  // Subtract the pedestal for normalization

      //running sum
      abs_run_sum(sample);
      sample = running_sum_;
     
      //running sum pedestal estimation + subtraction
      frugal_accum_update(sample, abs_pedestal_, abs_accum_);
      sample -= abs_pedestal_;
      //**no hit finding**

    }

    for (size_t i_t = 0; i_t < adcs.size(); ++i_t) {

      // if threshold < 0, the plane is not used to produce TPs
      if (threshold_ < 0)  continue;

      // get the sample
      int16_t sample = adcs[i_t];

      // update the pedestal estimate
      frugal_accum_update(sample, pedestal_, accum_);

      // pedestal subtract
      sample -= pedestal_;

      // apply the running sum
      abs_run_sum(sample);
      sample = running_sum_;

      // second pedestal estimation + subtraction
      frugal_accum_update(sample, abs_pedestal_, abs_accum_);
      sample -= abs_pedestal_;

      // check if we are over threshold
      bool is_over = sample > threshold_;
      if (is_over) {
        // we are over threshold, so need to update the hit charge and check for
        // peak time
        hit_charge_ += sample;

        // check if we're at the peak adc
        if (sample > hit_peak_adc_) {
          hit_peak_adc_ = (uint16_t)sample;
          hit_peak_time_ = hit_tover_;
        }

        // update time over threshold
        ++hit_tover_;
      }
      if (prev_was_over_ && !is_over) {
        // we've reached the end of the hit, so need to create a TP and write it
        // out

        this_tp.time_start = start_time + (i_t - hit_tover_) * this->ADC_SAMPLING_RATE_IN_DTS;
        this_tp.time_over_threshold = hit_tover_ * this->ADC_SAMPLING_RATE_IN_DTS;
        this_tp.time_peak = start_time + (i_t - hit_tover_ + hit_peak_time_) * this->ADC_SAMPLING_RATE_IN_DTS;
        this_tp.adc_integral = hit_charge_;
        this_tp.adc_peak = hit_peak_adc_;
        tps_out.push_back(this_tp);

        // now reset the hit variables
        hit_charge_ = 0;
        hit_tover_ = 0;
        hit_peak_adc_ = 0;
        hit_peak_time_ = 0;
      }
      prev_was_over_ = is_over;
    }
  }

private:
  // configuration parameters
  const int verbosity_;
  const int accum_limit_;

  // AbsRS-specific configuration parameters
  const int scale_factor_;
  const int r_value_;

  // Thresholds
  const int16_t threshold_tpg_plane0_;
  const int16_t threshold_tpg_plane1_;
  const int16_t threshold_tpg_plane2_;

  // variables for tracking state / hit-finding
  int16_t threshold_;
  int16_t pedestal_;
  int16_t accum_;

  int16_t abs_pedestal_;
  int16_t abs_accum_;
  int16_t running_sum_;
  // int16_t accum25_;
  //  int16_t accum75_;

  bool prev_was_over_;
  uint16_t hit_tover_;
  uint16_t hit_peak_time_;
  uint16_t hit_peak_adc_;
  uint32_t hit_charge_;
};

} // namespace dunetrigger

#endif
