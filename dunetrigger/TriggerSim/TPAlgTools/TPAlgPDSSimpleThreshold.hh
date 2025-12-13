#ifndef DUNETRIGGER_TRIGGERSIM_TPALGTPCSIMPLETHRESHOLD_hh
#define DUNETRIGGER_TRIGGERSIM_TPALGTPCSIMPLETHRESHOLD_hh

#include "fhiclcpp/ParameterSet.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/WireReadout.h"

#include "dunetrigger/TriggerSim/TPAlgTools/TPAlgTool.hh"
#include "dunetrigger/TriggerSim/Verbosity.hh"

#include <inttypes.h>
#include <limits>
#include <map>

namespace dunetrigger {

class TPAlgPDSSimpleThreshold : public TPAlgPDSTool {

public:
  explicit TPAlgPDSSimpleThreshold(fhicl::ParameterSet const &ps)
      : fVerbosity(ps.get<int>("verbosity", 0)),
        fPedWindowLength(ps.get<int>("pedestal_window_length", 20)),
        fThreshold(ps.get<int16_t>("threshold")) {}

  void
  initialize_channel_state(dunedaq::trgdataformats::channel_t const &channel,
                           std::vector<short> const &adcs) {
    if (fVerbosity >= Verbosity::kInfo) {
      std::cout << "Channel: " << channel << ", Threshold: " << fThreshold
                << std::endl;
    }

    // compute the pedestal from the first m_PedWindowLength samples
    size_t ped_length = std::min<size_t>(fPedWindowLength, adcs.size());
    fPedestal = std::accumulate(adcs.begin(), adcs.begin() + ped_length, 0.0);
    fPedestal /= ped_length;

    prev_was_over_ = false;
    hit_charge_ = 0;
    hit_tover_ = 0;
    hit_peak_adc_ = 0;
    hit_peak_time_ = 0;
  }

  void process_waveform(
      std::vector<short> const &adcs,
      dunedaq::trgdataformats::channel_t const channel,
      dunedaq::trgdataformats::detid_t const detid,
      dunedaq::trgdataformats::timestamp_t const start_time,
      std::vector<dunedaq::trgdataformats::TriggerPrimitive> &tps_out) {
    // setup a TP and initialize it with the common things for this
    // algorithm/channel
    dunedaq::trgdataformats::TriggerPrimitive this_tp;

    this_tp.channel = channel;
    this_tp.detid = detid;
    this_tp.type = dunedaq::trgdataformats::TriggerPrimitive::Type::kPDS;
    this_tp.algorithm =
        dunedaq::trgdataformats::TriggerPrimitive::Algorithm::kSimpleThreshold;
    this_tp.flag = 0;

    // for this channel, reinitialize the channel state variables
    initialize_channel_state(channel, adcs);

    for (size_t i_t = 0; i_t < adcs.size(); ++i_t) {

      // if threshold < 0, the plane is not used to produce TPs
      if (fThreshold < 0)
        continue;

      // get the sample
      int16_t sample = adcs[i_t] - fPedestal;

      // check if we are over threshold
      bool is_over = sample > fThreshold;
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

        this_tp.time_start =
            start_time + (i_t - hit_tover_) * this->ADC_SAMPLING_RATE_IN_DTS;
        this_tp.time_over_threshold =
            hit_tover_ * this->ADC_SAMPLING_RATE_IN_DTS;
        this_tp.time_peak = start_time + (i_t - hit_tover_ + hit_peak_time_) *
                                             this->ADC_SAMPLING_RATE_IN_DTS;
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
  const int fVerbosity;
  const int fPedWindowLength;
  const int16_t fThreshold;

  // variables for tracking state / hit-finding
  double fPedestal;

  bool prev_was_over_;
  uint16_t hit_tover_;
  uint16_t hit_peak_time_;
  uint16_t hit_peak_adc_;
  uint32_t hit_charge_;
};

} // namespace dunetrigger

#endif
