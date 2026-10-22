////////////////////////////////////////////////////////////////////////
// Class:       TriggerPrimitiveMakerTPC
// Plugin Type: producer (Unknown Unknown)
// File:        TriggerPrimitiveMakerTPC_module.cc
//
// Generated at Tue Nov 14 05:00:20 2023 by Wesley Ketchum using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

#include "dunetrigger/TriggerSim/TPAlgTools/TPAlgTPCTool.hh"

#include "detdataformats/DetID.hpp"
#include "detdataformats/trigger/TriggerPrimitive.hpp"
#include "lardataobj/RawData/RDTimeStamp.h"
#include "lardataobj/RawData/RawDigit.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/make_tool.h"
#include "canvas/Utilities/InputTag.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "canvas/Persistency/Common/FindOneP.h"
#include "dunetrigger/vendor/lardata/ArtDataHelper/GetManyByRegexTag.h"
#include "dunetrigger/TriggerSim/Verbosity.hh"
#include "dunetrigger/Utilities/ChannelOverlap.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "larcore/Geometry/WireReadout.h"

#include <algorithm>
#include <iostream>
#include <map>
#include <memory>
#include <unordered_map>
#include <utility>
#include <vector>

namespace dunetrigger {
class TriggerPrimitiveMakerTPC;
}

class dunetrigger::TriggerPrimitiveMakerTPC : public art::EDProducer {
public:
  explicit TriggerPrimitiveMakerTPC(fhicl::ParameterSet const &p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TriggerPrimitiveMakerTPC(TriggerPrimitiveMakerTPC const &) = delete;
  TriggerPrimitiveMakerTPC(TriggerPrimitiveMakerTPC &&) = delete;
  TriggerPrimitiveMakerTPC &
  operator=(TriggerPrimitiveMakerTPC const &) = delete;
  TriggerPrimitiveMakerTPC &operator=(TriggerPrimitiveMakerTPC &&) = delete;

  // Required functions.
  void produce(art::Event &e) override;

private:
  // Declare member data here.
  art::InputTag rawdigit_tag_;
  std::unique_ptr<TPAlgTPCTool> tpalg_;
  uint64_t default_timestamp_;
  int verbosity_;
};

dunetrigger::TriggerPrimitiveMakerTPC::TriggerPrimitiveMakerTPC(
    fhicl::ParameterSet const &p)
    : EDProducer{p} // ,
      ,
      rawdigit_tag_(p.get<art::InputTag>("rawdigit_tag")),
      tpalg_{art::make_tool<TPAlgTPCTool>(p.get<fhicl::ParameterSet>("tpalg"))},
      default_timestamp_(p.get<uint64_t>("default_timestamp", 0)),
      verbosity_(p.get<int>("verbosity", 0)) {
  // Call appropriate produces<>() functions here.
  // Call appropriate consumes<>() for any products to be retrieved by this
  // module.
  produces<std::vector<dunedaq::trgdataformats::TriggerPrimitive>>();
  consumesMany<std::vector<raw::RawDigit>>();
  consumesMany<art::Assns<raw::RDTimeStamp, raw::RawDigit>>();
}

void dunetrigger::TriggerPrimitiveMakerTPC::produce(art::Event &e) {
  // Implementation of required member function here.

  // make output collection for the TriggerPrimitive objects
  auto tp_col_ptr = std::make_unique<
      std::vector<dunedaq::trgdataformats::TriggerPrimitive>>();

  auto rawdigit_many =
      dunetrigger::vendor::lar::util::getManyByRegexTag<std::vector<raw::RawDigit>>(e, rawdigit_tag_);
  if (rawdigit_many.empty()) {
    throw cet::exception("TriggerPrimitiveMakerTPC")
        << "Found no std::vector<raw::RawDigit> collections matching rawdigit_tag \""
        << rawdigit_tag_.encode() << "\"";
  }

  // Collections sharing channels would silently produce duplicated TPs
  auto overlaps = channel_overlap::find(channel_overlap::from_handles(rawdigit_many));
  if (!overlaps.empty()) {
    channel_overlap::annotate(overlaps, art::ServiceHandle<geo::WireReadout>()->Get());
    throw cet::exception("TriggerPrimitiveMakerTPC")
        << "rawdigit_tag \"" << rawdigit_tag_.encode()
        << "\" matches raw::RawDigit collections sharing channels:\n"
        << channel_overlap::describe(overlaps)
        << "\nCheck that rawdigit_tag matches only one set of collections, e.g. by pinning the process name.";
  }

  for (auto const& rawdigit_handle : rawdigit_many) {

    // The timestamp associations are stored under the same tag as the
    // rawdigit collection they refer to.
    art::InputTag const rawdigit_tag = rawdigit_handle.provenance()->inputTag();

    if (verbosity_ >= Verbosity::kInfo)
      std::cout << "Processing " << rawdigit_tag.encode() << ": "
                << rawdigit_handle->size() << " raw::RawDigits" << std::endl;

    // try to get the associated timestamps to our rawdigit objects
    const art::FindOneP<raw::RDTimeStamp> rdtimestamp_per_rd(rawdigit_handle, e,
                                                            rawdigit_tag);
    // store a bool for whether it is valid or not to use inside the loop
    auto rd_assn_is_valid = rdtimestamp_per_rd.isValid();

    auto const& rawdigit_vec = *rawdigit_handle;

    uint64_t this_timestamp = default_timestamp_;
    for (size_t i_digit = 0; i_digit < rawdigit_vec.size(); ++i_digit) {
      auto const &digit = rawdigit_vec[i_digit];

      if (rd_assn_is_valid) {
        auto rdts = rdtimestamp_per_rd.at(i_digit);
        if (rdts)
          this_timestamp = rdts->GetTimeStamp();
      } else
        this_timestamp = default_timestamp_;

      tpalg_->process_waveform(
        digit.ADCs(), digit.Channel(),
        (uint16_t)(dunedaq::detdataformats::DetID::Subdetector::kHD_TPC),
        this_timestamp, *tp_col_ptr
      );

    }

  }

  e.put(std::move(tp_col_ptr));
}

DEFINE_ART_MODULE(dunetrigger::TriggerPrimitiveMakerTPC)
