////////////////////////////////////////////////////////////////////////
// Class:       TriggerCandidateMakerV5
// Plugin Type: producer (Unknown Unknown)
// File:        TriggerCandidateMakerV5_module.cc
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Persistency/Common/types.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "dunetrigger/TriggerSim/triggeralgs/v5/include/triggeralgs/TriggerCandidate.hpp"
#include "dunetrigger/TriggerSim/triggeralgs/v5/include/triggeralgs/TriggerCandidateFactory.hpp"

#include "detdataformats/trigger/TriggerActivityData2.hpp"
#include "detdataformats/trigger/TriggerCandidateData2.hpp"

#include "dunetrigger/TriggerSim/Verbosity.hh"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <memory>
#include <string>
#include <utility>

namespace dunedaq::trgdataformats2 {
bool operator==(const TriggerActivityData &ta1,
                const TriggerActivityData &ta2) {
  return (ta1.detid == ta2.detid && ta1.channel_peak == ta2.channel_peak &&
          ta1.time_peak == ta2.time_peak &&
          ta1.adc_integral == ta2.adc_integral);
}
bool operator<(const TriggerActivityData &ta1, const TriggerActivityData &ta2) {
  return std::tie(ta1.time_start, ta1.channel_start) <
         std::tie(ta2.time_start, ta2.channel_start);
}

} // namespace dunedaq::trgdataformats2

namespace dunetrigger {
class TriggerCandidateMakerV5;
using dunedaq::trgdataformats2::TriggerActivityData;
using dunedaq::trgdataformats2::TriggerCandidateData;
typedef std::pair<size_t, TriggerActivityData> TriggerActivityIdx;

} // namespace dunetrigger

class dunetrigger::TriggerCandidateMakerV5 : public art::EDProducer {
public:
  explicit TriggerCandidateMakerV5(fhicl::ParameterSet const &p);

  // Plugins should not be copied or assigned.
  TriggerCandidateMakerV5(TriggerCandidateMakerV5 const &) = delete;
  TriggerCandidateMakerV5(TriggerCandidateMakerV5 &&) = delete;
  TriggerCandidateMakerV5 &operator=(TriggerCandidateMakerV5 const &) = delete;
  TriggerCandidateMakerV5 &operator=(TriggerCandidateMakerV5 &&) = delete;

  void beginJob() override;
  void produce(art::Event &e) override;

private:
  art::InputTag ta_tag;
  std::string algname;
  fhicl::ParameterSet algconfig;

  std::shared_ptr<
      triggeralgs::AbstractFactory<triggeralgs::TriggerCandidateMaker>>
      alg_factory = triggeralgs::TriggerCandidateFactory::get_instance();
  std::unique_ptr<triggeralgs::TriggerCandidateMaker> alg;

  int verbosity;
};

dunetrigger::TriggerCandidateMakerV5::TriggerCandidateMakerV5(
    fhicl::ParameterSet const &p)
    : EDProducer{p}, ta_tag(p.get<art::InputTag>("ta_tag")),
      algname(p.get<std::string>("algorithm")),
      algconfig(p.get<fhicl::ParameterSet>("algconfig")),
      verbosity(p.get<int>("verbosity", 0)) {

  consumes<std::vector<TriggerActivityData>>(ta_tag);
  produces<std::vector<TriggerCandidateData>>();
  produces<art::Assns<TriggerCandidateData, TriggerActivityData>>();
}

void dunetrigger::TriggerCandidateMakerV5::beginJob() {
  // build alg using the factory
  alg = alg_factory->build_maker(algname);

  // parse the parameterset as json
  nlohmann::json alg_json;
  for (auto &k : algconfig.get_all_keys()) {
    try {
      alg_json[k] = algconfig.get<uint64_t>(k);
    } catch (const fhicl::exception &e) {
      try {
        // If false, try retrieving the parameter as a boolean
        alg_json[k] = algconfig.get<bool>(k);
      } catch (const fhicl::exception &e) {
        std::cerr << "Error: FHiCL parameter is neither an int nor a bool in "
                     "the FHiCL file. \n";
      }
    }
  }

  // and pass that on to the trigger algorithm
  alg->configure(alg_json);
}

void dunetrigger::TriggerCandidateMakerV5::produce(art::Event &e) {
  // get a handle on the TAs and dereference it
  auto ta_handle = e.getValidHandle<std::vector<TriggerActivityData>>(ta_tag);
  std::vector<TriggerActivityData> ta_vec = *ta_handle;

  // some unique pointers to new vectors for the data products
  auto tc_vec_ptr = std::make_unique<std::vector<TriggerCandidateData>>();
  auto ta_in_tc_assn_ptr =
      std::make_unique<art::Assns<TriggerCandidateData, TriggerActivityData>>();

  art::PtrMaker<TriggerCandidateData> tc_ptr_maker{e};

  // create a vector of inputs with the 'file' index of the TA
  std::vector<TriggerActivityIdx> input_tas; //(ta_vec.size());
  for (size_t i = 0; i < ta_vec.size(); ++i) {
    input_tas.push_back(TriggerActivityIdx(i, ta_vec.at(i)));
  }

  std::sort(input_tas.begin(), input_tas.end());

  // create a vector of online TCs for the online algorithm to store it's
  // outputs in
  std::vector<triggeralgs::TriggerCandidate> produced_tcs = {};

  // process the input TAs
  for (const auto &ta_indexed : input_tas) {
    TriggerActivityData ta_data = ta_indexed.second;
    triggeralgs::TriggerActivity curr_ta;
    static_cast<TriggerActivityData &>(curr_ta) = ta_data;
    (*alg)(curr_ta, produced_tcs);
  }

  // now we need to handle the associations
  for (auto out_tc : produced_tcs) {
    // make an art pointer to the tc vector
    auto const tcPtr = tc_ptr_maker(tc_vec_ptr->size());
    tc_vec_ptr->emplace_back(out_tc);

    // and now get the TAs out of the processed TC
    art::PtrVector<TriggerActivityData> tas_in_tc_ptr;
    for (const TriggerActivityData &in_ta : out_tc.inputs) {

      std::vector<TriggerActivityIdx>::iterator ta_it = std::find_if(
          input_tas.begin(), input_tas.end(),
          [&](const TriggerActivityIdx &curr_ta) { return curr_ta.second == in_ta; });
      // stop when there are no more matches
      while (ta_it != input_tas.end()) {
        tas_in_tc_ptr.push_back(
            art::Ptr<TriggerActivityData>(ta_handle, ta_it->first));
        // get an iterator pointing to the next match
        ta_it = std::find_if(++ta_it, input_tas.end(),
                             [&](const TriggerActivityIdx &curr_ta) {
                               return curr_ta.second == in_ta;
                             });
      }

      // print a debug message if we can't find and associated TAs
      if (tas_in_tc_ptr.empty() && verbosity >= Verbosity::kDebug) {
        std::cout << "No associated TAs found for TC!" << std::endl;
      }
    }
    // add the associations
    ta_in_tc_assn_ptr->addMany(tcPtr, tas_in_tc_ptr);
  }

  if (verbosity >= Verbosity::kInfo) {
    std::cout << "Created " << produced_tcs.size() << " TCs" << std::endl;
  }

  // move the produced things onto the event
  e.put(std::move(tc_vec_ptr));
  e.put(std::move(ta_in_tc_assn_ptr));
}

DEFINE_ART_MODULE(dunetrigger::TriggerCandidateMakerV5)
