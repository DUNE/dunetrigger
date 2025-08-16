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

  for ( auto mct : mctags ) {

    // Add this tag to the list of samples and get its index
    uint16_t mc_index = this->register_sample(mct);

    // Retrieve the truth vector for this tag
    auto mc_vec = e.getValidHandle<std::vector<simb::MCTruth>>(mct);

    // Loop over truths
    for( auto& mc : *mc_vec ) {
      // loop over mcparticles
      for( int k=0; k<mc.NParticles(); ++k) {
        // Map track ids to samples
        tk_to_sample.emplace(mc.GetParticle(k).TrackId(), mc_index);
      }
    }
  }

  return tk_to_sample;

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


  // std::map<uint16_t, const std::vector<simb::MCTruth>*> truths;
  // Loop over MC tags
  // std::cout << "xxxxxxxxxxxxxxxxxxx" << std::endl;
  // for ( auto t : mctags ) {

  //   // 1. Check if label exists in the index
  //   auto it = mc_label_to_index.find(t.label());
  //   uint16_t mc_index;
  //   if ( it != mc_label_to_index.end() ) {
  //     // 1.1 save the index
  //     mc_index = it->second;
  //   } else {
  //     // 1.2 or create a new one 
  //     mc_index = mc_label_to_index.size();

  //     // Populate the maps
  //     mc_label_to_index.emplace( t.label(), mc_index);
  //     mc_index_to_label.emplace( mc_index, t.label());
  //   }


    
  //   auto mc_vec = e.getValidHandle<std::vector<simb::MCTruth>>(t);

  //   std::cout << "label " << t.label() << " process " << t.process() << " instance " << t.instance() << " [" << t.encode() << "] : index = " << mc_index << std::endl;

  //   // Loop over the truth vector for this sample
  //   for( auto& mc : *mc_vec ) {
  //     // loop over particles
  //     for( int k=0; k<mc.NParticles(); ++k) {
  //       tk_to_sample_map[mc.GetParticle(k).TrackId()] = mc_index;
  //     }
  //   }

  // }
  // std::cout << "-------------------" << std::endl;


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
            std::cout << "WARN TID=0: simchan=" << sc.Channel() << ", ide=" << ide_i << ": t_id=" << ide.trackID << std::endl;
            continue;
        }

        if ( ide_i < 1000 ) {
          uint16_t sample_id = tk_to_sample_map[ide.trackID];
          std::cout << "simchan=" << sc.Channel() << ", ide=" << ide_i << ": t_id=" << ide.trackID << " ~ sample_id " << sample_id << " sample_name='" << mc_index_to_label[sample_id] << "'" << std::endl;
        }

      }
      ++ide_i;
      // fIDEsData.entries += tdcide.size();
    }
    ++sc_i;
  }

}

DEFINE_ART_MODULE(duneana::TriggerTricks)
