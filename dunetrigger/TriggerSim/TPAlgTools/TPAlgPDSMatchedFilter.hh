#ifndef DUNETRIGGER_TRIGGERSIM_TPALGPDSMATCHEDFILTER_hh
#define DUNETRIGGER_TRIGGERSIM_TPALGPDSMATCHEDFILTER_hh

#include "detdataformats/trigger/Types.hpp"
#include "fhiclcpp/ParameterSet.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/WireReadout.h"

#include "dunetrigger/TriggerSim/TPAlgTools/TPAlgTool.hh"
#include "dunetrigger/TriggerSim/Verbosity.hh"

#include <inttypes.h>

namespace dunetrigger {

class TPAlgPDSMatchedFilter : public TPAlgPDSTool {

public:
  explicit TPAlgPDSMatchedFilter(fhicl::ParameterSet const &ps)
      : fVerbosity(ps.get<int>("verbosity", 0)),
        fPedWindowLength(ps.get<int>("pedestal_window_length", 20)),
        fFixedPedestal(ps.get<int>("fixed_pedestal", -99999)),
        fThreshold(ps.get<float>("threshold")),
        fSubSampleFactor(ps.get<int>("subsample_factor", 3)),
        fTemplate(ps.get<std::vector<short>>("template")),
        fXCorrNormFactor(1.0 / std::accumulate(fTemplate.begin(),
                                               fTemplate.end(), 0.0,
                                               [](double sum, short val) {
                                                 return sum + val * val;
                                               })),
        fWindowDelay(ps.get<int>("window_delay", 0)),
        fWindowExtend(ps.get<int>("window_extend", 0)) {}

  void set_pedestal(std::vector<short> const &adcs) {
    if (fFixedPedestal != -99999) {
      fPedestal = fFixedPedestal;
    } else { // compute the pedestal from the first m_PedWindowLength samples
      size_t ped_length = std::min<size_t>(fPedWindowLength, adcs.size());
      fPedestal = std::accumulate(adcs.begin(), adcs.begin() + ped_length, 0.0);
      fPedestal /= ped_length;
    }
  }

  std::vector<double> correlate(std::vector<short> const &adcs);
  std::vector<std::pair<dunedaq::trgdataformats::timestamp_t,
                        dunedaq::trgdataformats::timestamp_t>>
  find_tp_windows(std::vector<double> const &xcorr);

  dunedaq::trgdataformats::TriggerPrimitive initalize_tp() const {
    dunedaq::trgdataformats::TriggerPrimitive result;
    result.channel = fCurrentChannel;
    result.detid = fCurrentDetID;
    result.type = dunedaq::trgdataformats::TriggerPrimitive::Type::kPDS;
    result.algorithm = dunedaq::trgdataformats::TriggerPrimitive::Algorithm::
        kSimpleThreshold; // FIXME
    result.flag = 0;
    return result;
  }

  void process_waveform(
      std::vector<short> const &adcs,
      dunedaq::trgdataformats::channel_t const channel,
      dunedaq::trgdataformats::detid_t const detid,
      dunedaq::trgdataformats::timestamp_t const start_time,
      std::vector<dunedaq::trgdataformats::TriggerPrimitive> &tps_out) override;

protected:
  // configuration parameters
  const int fVerbosity;
  const int fPedWindowLength;
  const int fFixedPedestal;
  const float fThreshold;
  const int fSubSampleFactor;
  const std::vector<short> fTemplate;
  const double fXCorrNormFactor;
  const int fWindowDelay;
  const int fWindowExtend;

  // variables for tracking state / hit-finding
  double fPedestal;
  dunedaq::trgdataformats::channel_t fCurrentChannel;
  dunedaq::trgdataformats::timestamp_t fCurrentDetID;
  dunedaq::trgdataformats::timestamp_t fCurrentStartTime;
};

} // namespace dunetrigger

#endif
