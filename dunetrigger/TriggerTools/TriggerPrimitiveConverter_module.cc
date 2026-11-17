////////////////////////////////////////////////////////////////////////
// Class:       TriggerPrimitiveConverter
// Plugin Type: producer
// File:        TriggerPrimitiveConverter_module.cc
//
// Converts dunedaq::trgdataformats::TriggerPrimitive (v1) collections into
// dunedaq::trgdataformats2::TriggerPrimitive (v2). The v2 format is packed,
// has narrower metadata fields, and expresses the peak and over-threshold
// times as ADC sample counts instead of DTS clock ticks. The v1 type and
// algorithm fields have no v2 counterpart and are dropped.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "cetlib_except/exception.h"

#include "detdataformats/trigger/TriggerPrimitive.hpp"
#include "detdataformats/trigger/TriggerPrimitive2.hpp"

#include <cstdint>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace dunetrigger {

class TriggerPrimitiveConverter : public art::EDProducer {
public:
  explicit TriggerPrimitiveConverter(fhicl::ParameterSet const& p);

  TriggerPrimitiveConverter(TriggerPrimitiveConverter const&) = delete;
  TriggerPrimitiveConverter(TriggerPrimitiveConverter&&) = delete;
  TriggerPrimitiveConverter& operator=(TriggerPrimitiveConverter const&) = delete;
  TriggerPrimitiveConverter& operator=(TriggerPrimitiveConverter&&) = delete;

  void produce(art::Event& e) override;
  void endJob() override;

private:
  using TriggerPrimitiveV1 = dunedaq::trgdataformats::TriggerPrimitive;
  using TriggerPrimitiveV2 = dunedaq::trgdataformats2::TriggerPrimitive;

  // Conversion counters, per event and summed over the job. The lossy counters
  // are valid v1 values that v2 cannot represent (warned about). invalid_fields
  // counts v1 sentinels mapped to the v2 sentinels (lossless, info only).
  struct ConversionStats {
    unsigned long flag_narrowed = 0;
    unsigned long detid_narrowed = 0;
    unsigned long channel_narrowed = 0;
    unsigned long sample_fields_clamped = 0;
    unsigned long sample_fields_rounded = 0;
    unsigned long peak_before_start = 0;
    unsigned long non_tpc = 0;
    unsigned long invalid_fields = 0;

    bool any_lossy() const
    {
      return flag_narrowed || detid_narrowed || channel_narrowed || sample_fields_clamped ||
             sample_fields_rounded || peak_before_start || non_tpc;
    }

    ConversionStats& operator+=(ConversionStats const& o)
    {
      flag_narrowed += o.flag_narrowed;
      detid_narrowed += o.detid_narrowed;
      channel_narrowed += o.channel_narrowed;
      sample_fields_clamped += o.sample_fields_clamped;
      sample_fields_rounded += o.sample_fields_rounded;
      peak_before_start += o.peak_before_start;
      non_tpc += o.non_tpc;
      invalid_fields += o.invalid_fields;
      return *this;
    }
  };

  TriggerPrimitiveV2 convert(TriggerPrimitiveV1 const& tp, ConversionStats& stats) const;
  // One "; "-separated clause per non-zero lossy counter.
  std::string describe_lossy(ConversionStats const& stats) const;
  uint16_t ticks_to_samples(dunedaq::trgdataformats::timestamp_t ticks,
                            uint16_t invalid_sentinel,
                            ConversionStats& stats) const;

  art::InputTag input_tag_;
  // DTS clock ticks per ADC sample: 32 for TPC TPs at 62.5 MHz DTS (ProtoDUNE-II,
  // DUNE FD). Other clocks or PDS sampling need a different value.
  dunedaq::trgdataformats::timestamp_t dts_ticks_per_sample_;
  // v2 has no type field, so non-TPC TPs would be indistinguishable from TPC
  // ones in the output. Reject them unless explicitly allowed.
  bool allow_non_tpc_;
  int verbosity_;

  // Job totals, reported in endJob().
  ConversionStats job_stats_;
  unsigned long job_events_ = 0;
  unsigned long job_events_lossy_ = 0;
  unsigned long job_tps_ = 0;
};

TriggerPrimitiveConverter::TriggerPrimitiveConverter(fhicl::ParameterSet const& p)
  : EDProducer{p}
  , input_tag_{p.get<art::InputTag>("input_tag")}
  , dts_ticks_per_sample_{p.get<dunedaq::trgdataformats::timestamp_t>("dts_ticks_per_sample", 32)}
  , allow_non_tpc_{p.get<bool>("allow_non_tpc", false)}
  , verbosity_{p.get<int>("verbosity", 0)}
{
  if (dts_ticks_per_sample_ == 0)
    throw cet::exception("TriggerPrimitiveConverter")
      << "dts_ticks_per_sample must be > 0 (got 0)\n";

  produces<std::vector<TriggerPrimitiveV2>>();
  consumes<std::vector<TriggerPrimitiveV1>>(input_tag_);
}

// Convert one valid DTS-tick interval into a count of ADC samples. Callers
// must filter out INVALID_TIMESTAMP inputs first. Intervals too long for the
// 16-bit v2 field become invalid_sentinel; intervals that are not a whole
// number of samples are rounded down. Both are counted in stats.
uint16_t TriggerPrimitiveConverter::ticks_to_samples(
  dunedaq::trgdataformats::timestamp_t ticks, uint16_t invalid_sentinel, ConversionStats& stats) const
{
  uint64_t const samples = ticks / dts_ticks_per_sample_;
  if (samples >= invalid_sentinel) {
    ++stats.sample_fields_clamped;
    return invalid_sentinel;
  }
  if (ticks % dts_ticks_per_sample_ != 0)
    ++stats.sample_fields_rounded;
  return static_cast<uint16_t>(samples);
}

TriggerPrimitiveConverter::TriggerPrimitiveV2
TriggerPrimitiveConverter::convert(TriggerPrimitiveV1 const& tp, ConversionStats& stats) const
{
  TriggerPrimitiveV2 out; // sets version = 2 and the invalid sentinels

  if (tp.type != TriggerPrimitiveV1::Type::kTPC)
    ++stats.non_tpc;

  // v1 Flags is 16 bits, v2 flag is 8: bits 8-15 are kUnassigned* today, so
  // keeping the low bits loses nothing in practice.
  out.flag = static_cast<uint8_t>(tp.flag & std::numeric_limits<uint8_t>::max());
  if (tp.flag > std::numeric_limits<uint8_t>::max())
    ++stats.flag_narrowed;

  if (tp.detid == dunedaq::trgdataformats::INVALID_DETID) {
    out.detid = dunedaq::trgdataformats2::INVALID_DETID;
    ++stats.invalid_fields;
  }
  else if (tp.detid == dunedaq::trgdataformats::WHOLE_DETECTOR)
    out.detid = dunedaq::trgdataformats2::WHOLE_DETECTOR;
  else if (tp.detid >= dunedaq::trgdataformats2::WHOLE_DETECTOR) {
    // Not representable in v2's uint8 detid: report as invalid rather than
    // silently truncating to a valid but wrong subdetector.
    out.detid = dunedaq::trgdataformats2::INVALID_DETID;
    ++stats.detid_narrowed;
  }
  else
    out.detid = static_cast<dunedaq::trgdataformats2::detid_t>(tp.detid);

  // v2 splits the 24-bit channel over three bytes (see Types2.hpp:
  // INVALID_TP_CHANNEL = 0xFFFFFF). Note the Channel() accessor in
  // TriggerPrimitive2.hpp mis-assembles the high byte (shifts by 12, not 16),
  // so the bytes are set directly here.
  dunedaq::trgdataformats2::channel_t channel = dunedaq::trgdataformats2::INVALID_TP_CHANNEL;
  if (tp.channel == dunedaq::trgdataformats::INVALID_CHANNEL)
    ++stats.invalid_fields;
  else if (tp.channel < 0 || tp.channel >= static_cast<dunedaq::trgdataformats::channel_t>(
                                             dunedaq::trgdataformats2::INVALID_TP_CHANNEL))
    ++stats.channel_narrowed;
  else
    channel = static_cast<dunedaq::trgdataformats2::channel_t>(tp.channel);

  // TODO: replace the three byte assignments below with
  // out.SetChannel(channel) once a dunedetdataformats release that includes
  // TriggerPrimitive::SetChannel (currently only on the develop branch) is
  // available in this area.
  out.channel_lowbyte = static_cast<uint8_t>(channel & 0xFF);
  out.channel_middlebyte = static_cast<uint8_t>((channel >> 8) & 0xFF);
  out.channel_highbyte = static_cast<uint8_t>((channel >> 16) & 0xFF);

  using dunedaq::trgdataformats::INVALID_TIMESTAMP;
  using dunedaq::trgdataformats2::INVALID_SAMPLES_OVER_THRESHOLD;
  using dunedaq::trgdataformats2::INVALID_SAMPLES_TO_PEAK;

  // Both formats use max(uint64) as INVALID_TIMESTAMP, so this is a plain copy.
  out.time_start = tp.time_start;

  // tp.time_over_threshold is a duration in DTS ticks.
  if (tp.time_over_threshold == INVALID_TIMESTAMP) {
    out.samples_over_threshold = INVALID_SAMPLES_OVER_THRESHOLD;
    ++stats.invalid_fields;
  }
  else
    out.samples_over_threshold =
      ticks_to_samples(tp.time_over_threshold, INVALID_SAMPLES_OVER_THRESHOLD, stats);

  // tp.time_peak is an absolute timestamp; v2 wants the offset from time_start
  // in samples. Check before subtracting: with an invalid time_start the
  // unsigned difference wraps to a small, plausible-looking value.
  if (tp.time_start == INVALID_TIMESTAMP || tp.time_peak == INVALID_TIMESTAMP) {
    out.samples_to_peak = INVALID_SAMPLES_TO_PEAK;
    ++stats.invalid_fields;
  }
  else if (tp.time_peak < tp.time_start) {
    out.samples_to_peak = INVALID_SAMPLES_TO_PEAK;
    ++stats.peak_before_start;
  }
  else
    out.samples_to_peak =
      ticks_to_samples(tp.time_peak - tp.time_start, INVALID_SAMPLES_TO_PEAK, stats);

  out.adc_integral = tp.adc_integral;
  out.adc_peak = tp.adc_peak;

  return out;
}

void TriggerPrimitiveConverter::produce(art::Event& e)
{
  auto const& tps_v1 = *e.getValidHandle<std::vector<TriggerPrimitiveV1>>(input_tag_);

  auto tps_v2 = std::make_unique<std::vector<TriggerPrimitiveV2>>();
  tps_v2->reserve(tps_v1.size());

  ConversionStats stats;
  for (auto const& tp : tps_v1)
    tps_v2->push_back(convert(tp, stats));

  if (stats.non_tpc && !allow_non_tpc_)
    throw cet::exception("TriggerPrimitiveConverter")
      << input_tag_.encode() << ": " << stats.non_tpc << " of " << tps_v1.size()
      << " TP(s) are not Type::kTPC. v2 has no type field, so they would be"
      << " indistinguishable from TPC TPs. Set allow_non_tpc: true to convert anyway.\n";

  ++job_events_;
  job_tps_ += tps_v1.size();
  job_stats_ += stats;

  if (stats.any_lossy()) {
    ++job_events_lossy_;
    mf::LogWarning("TriggerPrimitiveConverter")
      << "Narrowed or clamped fields converting " << input_tag_.encode() << ":"
      << describe_lossy(stats);
  }

  if (verbosity_ >= 1)
    mf::LogInfo("TriggerPrimitiveConverter")
      << input_tag_.encode() << ": converted " << tps_v1.size() << " TP(s) to v2, "
      << stats.invalid_fields << " invalid input field(s) mapped to v2 sentinels";

  e.put(std::move(tps_v2));
}

std::string TriggerPrimitiveConverter::describe_lossy(ConversionStats const& stats) const
{
  std::ostringstream msg;
  if (stats.flag_narrowed)
    msg << " flag>255 in " << stats.flag_narrowed << " TP(s), kept low 8 bits;";
  if (stats.detid_narrowed)
    msg << " unrepresentable detid in " << stats.detid_narrowed << " TP(s), set to INVALID_DETID;";
  if (stats.channel_narrowed)
    msg << " invalid/out-of-range channel in " << stats.channel_narrowed
        << " TP(s), set to INVALID_TP_CHANNEL;";
  if (stats.sample_fields_clamped)
    msg << " over-long time interval in " << stats.sample_fields_clamped
        << " field(s), set to INVALID_SAMPLES_*;";
  if (stats.sample_fields_rounded)
    msg << " interval not a multiple of " << dts_ticks_per_sample_ << " DTS ticks in "
        << stats.sample_fields_rounded << " field(s), rounded down to whole samples;";
  if (stats.non_tpc)
    msg << " non-TPC type in " << stats.non_tpc << " TP(s), converted (allow_non_tpc) but"
        << " the type is lost;";
  if (stats.peak_before_start)
    msg << " time_peak before time_start in " << stats.peak_before_start
        << " TP(s), set to INVALID_SAMPLES_TO_PEAK;";
  return msg.str();
}

void TriggerPrimitiveConverter::endJob()
{
  std::ostringstream msg;
  msg << "Job summary for " << input_tag_.encode() << ": converted " << job_tps_
      << " TP(s) in " << job_events_ << " event(s); " << job_stats_.invalid_fields
      << " invalid input field(s) mapped to v2 sentinels.";
  if (job_stats_.any_lossy()) {
    msg << " Lossy conversions in " << job_events_lossy_ << " event(s):"
        << describe_lossy(job_stats_);
    mf::LogWarning("TriggerPrimitiveConverter") << msg.str();
  }
  else {
    msg << " No lossy conversions.";
    mf::LogInfo("TriggerPrimitiveConverter") << msg.str();
  }
}

} // namespace dunetrigger

DEFINE_ART_MODULE(dunetrigger::TriggerPrimitiveConverter)
