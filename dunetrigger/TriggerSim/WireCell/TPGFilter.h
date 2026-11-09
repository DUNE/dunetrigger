/** wclsTPGFilter

    A Wire-Cell frame filter that runs a dunetrigger TPAlgTPC tool on every
    channel of the (raw, digitized) frame of one anode/CRM and hands the
    resulting trigger primitives to a wclsTPGSaver.

    Place one instance per anode directly after the Digitizer, before any
    noise filtering or signal processing.

    Each channel's waveform is built exactly as wclsFrameSaver builds a
    raw::RawDigit (first trace of the channel, dense from tick 0, samples
    truncated to short, non-finite samples set to 0), so the TPs are
    identical to running TriggerPrimitiveMakerTPC on those RawDigits.

    Configuration:
      anode              "type:name" of the IAnodePlane (sets channel -> plane)
      saver              "type:name" of the wclsTPGSaver (default "wclsTPGSaver")
      tpalg              object: a TPAlgTPC* tool configuration, as it would be
                         written in FHiCL (must contain "tool_type")
      detid              TP detid (default 3 = kHD_TPC, same as
                         TriggerPrimitiveMakerTPC)
      default_timestamp  timestamp (DTS ticks) of sample 0 (default 0, same
                         as TriggerPrimitiveMakerTPC without RDTimeStamps)
      nticks             if > 0, force this waveform length (as wclsFrameSaver)
      scale              multiply samples before truncation (default 1.0)
      pass_through       true: output the input frame unchanged (needed if
                         signal processing follows). false: output an empty
                         frame with the same ident/time to free memory early.
*/

#ifndef DUNETRIGGER_TRIGGERSIM_WIRECELL_TPGFILTER_H
#define DUNETRIGGER_TRIGGERSIM_WIRECELL_TPGFILTER_H

#include "WireCellAux/Logger.h"
#include "WireCellIface/IAnodePlane.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IFrameFilter.h"

#include "dunetrigger/TriggerSim/TPAlgTools/TPAlgTPCTool.hh"

#include <cstdint>
#include <memory>
#include <mutex>
#include <unordered_map>

namespace dunetrigger {

  class TPGSaver;

  class TPGFilter : public WireCell::Aux::Logger,
                    public WireCell::IFrameFilter,
                    public WireCell::IConfigurable {
  public:
    TPGFilter();
    virtual ~TPGFilter();

    // IFrameFilter
    bool operator()(const input_pointer& in, output_pointer& out) override;

    // IConfigurable
    void configure(const WireCell::Configuration& config) override;
    WireCell::Configuration default_configuration() const override;

  private:
    WireCell::IAnodePlane::pointer m_anode;
    std::shared_ptr<TPGSaver> m_saver;
    std::unique_ptr<TPAlgTPCTool> m_tpalg;
    std::unordered_map<int, unsigned int> m_plane; // channel ident -> plane index

    int m_detid{3};
    uint64_t m_default_timestamp{0};
    int m_nticks{0};
    double m_scale{1.0};
    bool m_pass_through{true};

    std::mutex m_mutex; // the TP tool keeps per-channel state
    size_t m_count{0};
  };

}

#endif
