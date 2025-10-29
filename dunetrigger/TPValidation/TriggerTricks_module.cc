////////////////////////////////////////////////////////////////////////
// Class:       TriggerTricks
// Plugin Type: analyzer (Unknown Unknown)
// File:        TriggerTricks_module.cc
//
// Generated at Fri Aug 15 14:35:52 2025 by Alessandro Thea,1 1-006,+41227666345, using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"

#include "lardataobj/Simulation/SimChannel.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

namespace duneana {
  class TriggerTricks;
}


class duneana::TriggerTricks : public art::EDAnalyzer {
public:
  explicit TriggerTricks(fhicl::ParameterSet const& p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TriggerTricks(TriggerTricks const&) = delete;
  TriggerTricks(TriggerTricks&&) = delete;
  TriggerTricks& operator=(TriggerTricks const&) = delete;
  TriggerTricks& operator=(TriggerTricks&&) = delete;

  // Required functions.
  void analyze(art::Event const& e) override;

  // Update the sample maps with 
  uint16_t register_sample(const art::InputTag& mctag);
  
  // Build the track to sample map 
  std::map<int, uint16_t> make_tk_sample_map(art::Event const& e);
  std::map<art::Ptr<simb::MCTruth>, uint16_t> make_mctruth_to_sample_map(art::Event const& e);

  // Particle Inventory Service
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;

  // MC Sample to index map
  std::map<std::string, uint16_t> mc_label_to_index;
  // MC index to sample map
  std::map<uint16_t,std::string> mc_index_to_label;


private:

  // Declare member data here.

};


duneana::TriggerTricks::TriggerTricks(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}


uint16_t 
duneana::TriggerTricks::register_sample(const art::InputTag& mc_tag) {

  // Update the mc samples maps
    // 1. Check if label exists in the index
    auto it = mc_label_to_index.find(mc_tag.label());
    uint16_t mc_index;
    if ( it != mc_label_to_index.end() ) {
      // 1.1 This samp
      mc_index = it->second;
    } else {
      // 1.2 or create a new one 
      mc_index = mc_label_to_index.size();

      // Populate the maps
      mc_label_to_index.emplace( mc_tag.label(), mc_index);
      mc_index_to_label.emplace( mc_index, mc_tag.label());
    }

    return mc_index;
}

std::map<int, uint16_t>
duneana::TriggerTricks::make_tk_sample_map(art::Event const& e) {

  auto mctags = e.getInputTags<std::vector<simb::MCTruth>>();

  std::map<int, uint16_t> tk_to_sample;

  std::cout << "--- MAKING MAP ---" << std::endl;


  std::map<art::Ptr<simb::MCTruth>, uint16_t> mct_map;

  for ( auto mct : mctags ) {

    // Add this tag to the list of samples and get its index
    uint16_t mc_index = this->register_sample(mct);

    // Retrieve the truth vector for this tag
    auto mc_vec = e.getValidHandle<std::vector<simb::MCTruth>>(mct);
    std::cout << "--- Truth ---" << mct << std::endl;

    // Loop over truths
    for( auto& mc : *mc_vec ) {
      // loop over mcparticles
      for( int k=0; k<mc.NParticles(); ++k) {
        // Map track ids to samples
        tk_to_sample.emplace(mc.GetParticle(k).TrackId(), mc_index);
        std::cout << "+++ track id " << mc.GetParticle(k).TrackId() << " -> " << mc_index << std::endl;
      }
    }
  }
  std::cout << "--- MAKING MAP ---" << std::endl;

  return tk_to_sample;

}



std::map<art::Ptr<simb::MCTruth>, uint16_t>
duneana::TriggerTricks::make_mctruth_to_sample_map(art::Event const& e) {

  auto mctags = e.getInputTags<std::vector<simb::MCTruth>>();

  std::cout << "--- MAKING TRUTH MAP ---" << std::endl;

  std::map<art::Ptr<simb::MCTruth>, uint16_t> mct_to_sample_map;

  for ( auto mct : mctags ) {

    // Add this tag to the list of samples and get its index
    uint16_t mc_index = this->register_sample(mct);

    // Retrieve the truth vector for this tag
    auto mc_vec = e.getValidHandle<std::vector<simb::MCTruth>>(mct);
    std::cout << "--- Truth ---" << mct << std::endl;
    std::cout << "    ProdId : " << mc_vec.id() << std::endl;

    for(size_t i(0); i<mc_vec->size(); ++i) {
      art::Ptr<simb::MCTruth> ptr(mc_vec, i);
      mct_to_sample_map[ptr] = mc_index;
      std::cout << ptr << " ---- " << mc_index << std::endl;
    }

  }
  std::cout << "--- MAKING TRUTH MAP ---" << std::endl;

  return mct_to_sample_map;

}


void 
duneana::TriggerTricks::analyze(art::Event const& e)
{
  // Implementation of required member function here.


  std::cout << "Analyzing!" << std::endl;
  std::cout << "Number of particles in the inventory: " << pi_serv->ParticleList().size() << std::endl;
  std::cout << "Number of mcthruths inventory: " << pi_serv->MCTruthVector_Ps().size() << std::endl;

  auto simchannels = e.getValidHandle<std::vector<sim::SimChannel>>("tpcrawdecoder:simpleSC");

  std::map<int, uint16_t> tk_to_sample_map = this->make_tk_sample_map(e);
  auto mctp_to_sample = make_mctruth_to_sample_map(e);

  for( auto& [tk_id, mc_idx] : tk_to_sample_map ) {
    std::cout << "--- track id " << tk_id << " -> " << mc_idx << std::endl;
  }

  for( auto& [label, mc_index] : mc_label_to_index ) {
    std::cout << "label = " << label << ", index = " << mc_index << std::endl;
  }

  size_t sc_i(0), ide_i(0);

  // Process the SimChannels just once
  for (const sim::SimChannel &sc : *simchannels) {
    auto tdcidemap = sc.TDCIDEMap();

    // Sum charge for each plane
    for (const auto& [t, tdcide] : tdcidemap) {
      for (const auto& ide : tdcide) {
        if ( ide.trackID == 0) {
            std::cout << "WARN TID=0: simchan=" << sc.Channel() << ", ide=" << ide_i << ": t_id=" << ide.trackID << " orig tk_id: " <<  ide.origTrackID << std::endl;
            continue;
        }

        if ( ide_i < 1000 ) {
          // uint16_t sample_id = tk_to_sample_map.at(ide.trackID);

          const art::Ptr< simb::MCTruth >& ptr = pi_serv->TrackIdToMCTruth_P (ide.trackID);

          std::cout << "simchan=" << sc.Channel() << ", ide=" << ide_i << ", t_id=" << ide.trackID << " ~~~ prd id: " << ptr.id() << " ~~~ mct idx: " << mctp_to_sample.at(ptr) << std::endl;
        }

      }
      ++ide_i;
      // fIDEsData.entries += tdcide.size();
    }
    ++sc_i;
  }

}

DEFINE_ART_MODULE(duneana::TriggerTricks)
