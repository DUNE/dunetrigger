////////////////////////////////////////////////////////////////////////
// Class:       TriggerPrimitiveMakerPDS
// Plugin Type: producer (Unknown Unknown)
// File:        TriggerPrimitiveMakerPDS_module.cc
//
// Generated at Fri Dec 12 15:37:01 2025 by jierans using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Utilities/InputTag.h"
#include "detdataformats/trigger/Types.hpp"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "detdataformats/DetID.hpp"
#include "detdataformats/trigger/TriggerPrimitive.hpp"
#include "dunetrigger/TriggerSim/TPAlgTools/TPAlgTool.hh"
#include "dunetrigger/TriggerSim/Verbosity.hh"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RawData/OpDetWaveform.h"

#include <memory>

namespace dunetrigger {
class TriggerPrimitiveMakerPDS;
}

class dunetrigger::TriggerPrimitiveMakerPDS : public art::EDProducer {
public:
  explicit TriggerPrimitiveMakerPDS(fhicl::ParameterSet const &p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TriggerPrimitiveMakerPDS(TriggerPrimitiveMakerPDS const &) = delete;
  TriggerPrimitiveMakerPDS(TriggerPrimitiveMakerPDS &&) = delete;
  TriggerPrimitiveMakerPDS &
  operator=(TriggerPrimitiveMakerPDS const &) = delete;
  TriggerPrimitiveMakerPDS &operator=(TriggerPrimitiveMakerPDS &&) = delete;

  // Required functions.
  void produce(art::Event &e) override;

private:
  // Declare member data here.
  art::InputTag fOptDetWaveformTag;
  std::unique_ptr<TPAlgPDSTool> fTPAlg;
  int fVerbosity;
};

dunetrigger::TriggerPrimitiveMakerPDS::TriggerPrimitiveMakerPDS(
    fhicl::ParameterSet const &p)
    : EDProducer{p},
      fOptDetWaveformTag(p.get<art::InputTag>("opdetwaveform_tag")),
      fTPAlg{art::make_tool<TPAlgPDSTool>(p.get<fhicl::ParameterSet>("tpalg"))},
      fVerbosity(p.get<int>("verbosity", Verbosity::kInfo)) {
  // Call appropriate produces<>() functions here.
  // Call appropriate consumes<>() for any products to be retrieved by this
  // module.
  produces<std::vector<dunedaq::trgdataformats::TriggerPrimitive>>();
  consumes<std::vector<raw::OpDetWaveform>>(fOptDetWaveformTag);
}

void dunetrigger::TriggerPrimitiveMakerPDS::produce(art::Event &e) {
  auto const detectorClockData =
      art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
  auto tp_col_ptr = std::make_unique<
      std::vector<dunedaq::trgdataformats::TriggerPrimitive>>();
  auto opdetwaveform_handle =
      e.getHandle<std::vector<raw::OpDetWaveform>>(fOptDetWaveformTag);
  if (fVerbosity >= Verbosity::kInfo) {
    std::cout << "TriggerPrimitiveMakerPDS: found "
              << opdetwaveform_handle->size() << " OpDetWaveform objects."
              << std::endl;
  }
  for (auto wfm : *opdetwaveform_handle) {
    int timestamp_in_ticks =
        wfm.TimeStamp() / detectorClockData.OpticalClock().TickPeriod();
    // set start_time to zero here to avoid a negative timiestamp being
    // converted to a uint
    fTPAlg->process_waveform(
        wfm.Waveform(),
        static_cast<dunedaq::trgdataformats::channel_t>(wfm.ChannelNumber()),
        static_cast<dunedaq::trgdataformats::detid_t>(
            dunedaq::detdataformats::DetID::Subdetector::kHD_PDS),
        timestamp_in_ticks, *tp_col_ptr);
  }
  if (fVerbosity >= Verbosity::kInfo) {
    std::cout << "TriggerPrimitiveMakerPDS: Generated " << tp_col_ptr->size()
              << " TPs." << std::endl;
  }
  e.put(std::move(tp_col_ptr));
}

DEFINE_ART_MODULE(dunetrigger::TriggerPrimitiveMakerPDS)
