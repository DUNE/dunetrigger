#include "art/Utilities/ToolMacros.h"

#include "dunetrigger/TriggerSim/TPAlgTools/TPAlgPDSMatchedFilter.hh"

namespace dunetrigger {
using namespace dunedaq::trgdataformats;

std::vector<double>
TPAlgPDSMatchedFilter::correlate(std::vector<double> const &adcs) {
  size_t n_samples = adcs.size() / fSubSampleFactor;
  size_t n_taps = fTemplate.size();
  size_t n_out = n_samples - n_taps + 1;
  std::vector<double> xcorr(n_out, 0.0);

  for (size_t i = 0; i < n_out; i++) {
    double sum = 0.0;
    for (size_t j = 0; j < n_taps; ++j) {
      double adc_val = adcs.at(fSubSampleFactor * (i + j));
      sum += adc_val * static_cast<double>(fTemplate[j]);
    }
    xcorr[i] = sum * fXCorrNormFactor;
  }
  return xcorr;
}

std::vector<TPWindow>
TPAlgPDSMatchedFilter::find_tp_windows(std::vector<double> const &xcorr) {
  std::vector<TPWindow> tp_windows;
  bool window_open = false;
  int extend_counter = 0;
  for (size_t i = 0; i < xcorr.size(); ++i) {
    timestamp_t current_time = i * fSubSampleFactor + fWindowDelay;
    if (xcorr[i] >= fThreshold) {
      if (!window_open) {
        // start a new window
        tp_windows.emplace_back(current_time, current_time, 0.0);
      }
      window_open = true;
      std::get<1>(tp_windows.back()) = current_time;
      std::get<2>(tp_windows.back()) =
          std::max(std::get<2>(tp_windows.back()), xcorr[i]);
      extend_counter = 0;
    } else {
      if (window_open) {
        if (extend_counter < fWindowExtend) {
          // extend the current window
          std::get<1>(tp_windows.back()) = current_time;
          std::get<2>(tp_windows.back()) =
              std::max(std::get<2>(tp_windows.back()), xcorr[i]);
          extend_counter += fSubSampleFactor;
        } else {
          // close the current window
          window_open = false;
        }
      }
    }
  }
  return tp_windows;
}

void TPAlgPDSMatchedFilter::process_waveform(
    std::vector<short> const &adcs, channel_t const channel,
    detid_t const detid, timestamp_t const start_time,
    std::vector<TriggerPrimitive> &tps_out) {

  fCurrentChannel = channel;
  fCurrentDetID = detid;
  fCurrentStartTime = start_time;
  if (fVerbosity >= Verbosity::kInfo) {
    std::cout << "Channel: " << channel << ", Threshold: " << fThreshold
              << std::endl;
  }
  set_pedestal(adcs);
  std::vector<double> ped_subtracted(adcs.size(), 0.0);
  for (size_t i = 0; i < adcs.size(); ++i) {
    ped_subtracted[i] = static_cast<double>(adcs.at(i)) - fPedestal;
  }
  std::vector<double> xcorr = correlate(ped_subtracted);
  std::vector<TPWindow> tp_windows = find_tp_windows(xcorr);
  for (auto const &[start, end, npes] : tp_windows) {
    TriggerPrimitive this_tp = initalize_tp();
    this_tp.time_start = start_time + start;
    this_tp.time_over_threshold = (end - start + 1);
    this_tp.adc_peak = 0;
    double integral = 0.0;
    for (timestamp_t t = start; t <= end; ++t) {
      integral += ped_subtracted.at(t);
      if (ped_subtracted.at(t) > this_tp.adc_peak) {
        this_tp.adc_peak = ped_subtracted.at(t);
        this_tp.time_peak = start_time + t;
      }
    }
    this_tp.adc_integral = static_cast<uint32_t>(integral);
    this_tp.time_start *= ADC_SAMPLING_RATE_IN_DTS;
    this_tp.time_peak *= ADC_SAMPLING_RATE_IN_DTS;
    this_tp.time_over_threshold *= ADC_SAMPLING_RATE_IN_DTS;
    this_tp.adc_peak = static_cast<uint16_t>(npes * 10.0);
    tps_out.push_back(this_tp);
  }
}

} // namespace dunetrigger

DEFINE_ART_CLASS_TOOL(dunetrigger::TPAlgPDSMatchedFilter)
