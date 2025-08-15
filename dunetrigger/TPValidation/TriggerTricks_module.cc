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
#include "messagefacility/MessageLogger/MessageLogger.h"

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

private:

  // Declare member data here.

};


duneana::TriggerTricks::TriggerTricks(fhicl::ParameterSet const& p)
  : EDAnalyzer{p}  // ,
  // More initializers here.
{
  // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void duneana::TriggerTricks::analyze(art::Event const& e)
{
  // Implementation of required member function here.

  std::cout << "Analyzing!" << std::endl;
}

DEFINE_ART_MODULE(duneana::TriggerTricks)
