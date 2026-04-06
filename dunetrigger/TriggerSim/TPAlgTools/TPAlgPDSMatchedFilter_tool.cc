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

std::vector<MFWindow>
TPAlgPDSMatchedFilter::make_mf_windows(std::vector<double> const &xcorr) {
  std::vector<MFWindow> mf_windows;
  bool window_open = false;
  int extend_counter = 0;
  for (size_t i = 0; i < xcorr.size(); ++i) {
    timestamp_t current_time = i * fSubSampleFactor + fWindowDelay;
    bool above_threshold = (xcorr[i] >= fThreshold);
    if (!window_open && above_threshold) {
      // start a new window
      mf_windows.emplace_back(current_time);
      window_open = true;
    }
    if (window_open && (above_threshold || extend_counter < fWindowExtend)) {
      // extend the current window
      MFWindow &last_window = mf_windows.back();
      last_window.end_time = current_time;
      last_window.integral += xcorr[i];
      if (xcorr[i] > last_window.peak_value) {
        last_window.peak_time = current_time;
        last_window.peak_value = xcorr[i];
      }
      if (above_threshold) {
        extend_counter = 0;
      } else {
        extend_counter += fSubSampleFactor;
      }
    }
    if (window_open && !above_threshold && extend_counter >= fWindowExtend) {
      // close existing window
      window_open = false;
    }
  }
  return mf_windows;
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
  std::vector<MFWindow> tp_windows = make_mf_windows(xcorr);
  for (auto const &mf_window : tp_windows) {
    TriggerPrimitive this_tp = initalize_tp();
    this_tp.time_start = start_time + mf_window.start_time;
    this_tp.time_over_threshold =
        (mf_window.end_time - mf_window.start_time + 1);
    this_tp.time_start *= ADC_SAMPLING_RATE_IN_DTS;
    this_tp.time_over_threshold *= ADC_SAMPLING_RATE_IN_DTS;

    // --- Using raw waveform for peak and integral calculation ---
    this_tp.adc_peak = 0;
    double integral = 0.0;
    for (timestamp_t t = mf_window.start_time; t <= mf_window.end_time; ++t)
    {
      integral += ped_subtracted.at(t);
      if (ped_subtracted.at(t) > this_tp.adc_peak) {
        this_tp.adc_peak = ped_subtracted.at(t);
        this_tp.time_peak = start_time + t;
      }
    }
    this_tp.adc_integral = static_cast<uint32_t>(integral);

    // --- Using matched filter output for peak and integral calculation ---
    // this_tp.adc_integral = static_cast<uint32_t>(mf_window.integral * 10.0);
    this_tp.time_peak = mf_window.peak_time;
    this_tp.adc_peak = static_cast<uint16_t>(mf_window.peak_value * 10.0);
    // ---
    this_tp.time_peak *= ADC_SAMPLING_RATE_IN_DTS;
    tps_out.push_back(this_tp);
  }
}

} // namespace dunetrigger

DEFINE_ART_CLASS_TOOL(dunetrigger::TPAlgPDSMatchedFilter)
