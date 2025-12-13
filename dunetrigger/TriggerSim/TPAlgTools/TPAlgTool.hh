#ifndef DUNETRIGGER_TRIGGERSIM_TPALGTPCTOOL_hh
#define DUNETRIGGER_TRIGGERSIM_TPALGTPCTOOL_hh

#include "detdataformats/trigger/TriggerPrimitive.hpp"

#include <cstdint>
#include <vector>

namespace dunetrigger {

class TPAlgTool {
public:
  virtual ~TPAlgTool() noexcept = default;
  // take in a waveform, add trigger primitives to it
  virtual void process_waveform(
      std::vector<short> const &adcs,
      dunedaq::trgdataformats::channel_t const channel,
      dunedaq::trgdataformats::detid_t const detid,
      dunedaq::trgdataformats::timestamp_t const start_time,
      std::vector<dunedaq::trgdataformats::TriggerPrimitive> &tps_out) = 0;
};

class TPAlgTPCTool : public TPAlgTool {

public:
  static const unsigned int ADC_SAMPLING_RATE_IN_DTS =
      32; // 32 DTS time ticks between adc samples

  inline int16_t avx2_divide(const int16_t &a, const int16_t &b) {
    int16_t vb = (1 << 15) / b;  //  1 / b * 2^15
    int32_t mulhrs = a * vb;     //  a / b * 2^15
    mulhrs = (mulhrs >> 14) + 1; // (a / b * 2^15) * 2^-14 + 1 ~ a / b * 2 + 1
    mulhrs = mulhrs >> 1;        //~ a / b. The +1 causes unorthodox rounding.
    return (int16_t)(mulhrs);
  }
};

class TPAlgPDSTool : public TPAlgTool {
public:
  static const unsigned int ADC_SAMPLING_RATE_IN_DTS =
      1; // 1 DTS time tick between adc samples
};
} // namespace dunetrigger

#endif
