////////////////////////////////////////////////////////////////////////
// Class:       HydraFilter
// Plugin Type: filter
// File:        HydraFilter_module.cc
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/WireReadout.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"

#include <memory>

namespace dunetrigger {
class HydraFilter;
}

class dunetrigger::HydraFilter : public art::EDFilter {
public:
  explicit HydraFilter(fhicl::ParameterSet const& p);

  HydraFilter(HydraFilter const&) = delete;
  HydraFilter(HydraFilter&&) = delete;
  HydraFilter& operator=(HydraFilter const&) = delete;
  HydraFilter& operator=(HydraFilter&&) = delete;

  bool filter(art::Event& e) override;

private:
  short unsigned int tpcset_id;
};

dunetrigger::HydraFilter::HydraFilter(fhicl::ParameterSet const& p)
    : EDFilter{p},
    tpcset_id{p.get<short unsigned int>("tpcset_id")}
{
  // get a service handle for geometry
  geo::WireReadoutGeom const *wrgeom = &art::ServiceHandle<geo::WireReadout>()->Get();
  std::cout << "maxtpcsets = " <<  wrgeom->MaxTPCsets() << std::endl;

  // Check that tpcset_id exists
  auto tpc_vec = wrgeom->TPCsetToTPCs({0, tpcset_id});
  std::cout << "Num TPCs in TPCSet " << tpcset_id << " : " << tpc_vec.size() << std::endl;

}

bool dunetrigger::HydraFilter::filter(art::Event& e)
{

  auto const &geo = art::ServiceHandle<geo::Geometry>();
  std::cout << "Detector name: " << geo->DetectorName() << std::endl;


  // get a service handle for geometry
  geo::WireReadoutGeom const *wrgeom = &art::ServiceHandle<geo::WireReadout>()->Get();
  // std::cout << "maxtpcsets = " <<  wrgeom->MaxTPCsets() << std::endl;
  auto tpcid_vec = wrgeom->TPCsetToTPCs({0, tpcset_id});

  for ( auto tpcid : tpcid_vec ) {
    geo->TPC(tpcid);
  }

  


  // std::cout << "Num TPCs in TPCSet 0 : " << tpcid_vec.size() << std::endl;


  std::string m_inputTag = "IonAndScint";
  art::Handle<std::vector<sim::SimEnergyDeposit>> sedvh;
  bool okay = e.getByLabel(m_inputTag, sedvh);

  if ( !okay ) return false;

  // std::cout << okay << std::endl;

  bool has_deposits_in_tpcset = false;

  for( auto& tpcid: tpcid_vec) {
    auto& tpcgeo = geo->TPC(tpcid);
    // std::cout << "Testing TPC " << tpcid << std::endl;
  // for( unsigned int tpcid(0); tpcid < 112; ++tpcid) {
    // auto& tpcgeo = geo->TPC( geo::TPCID{0, tpcid});

    uint32_t num_deposits_in_tpc = false;
    for( auto& sed : *sedvh) {
      num_deposits_in_tpc += (
        sed.NumElectrons() > 0 &&
        (
          tpcgeo.ActiveBoundingBox().ContainsPosition(sed.Start()) ||
          tpcgeo.ActiveBoundingBox().ContainsPosition(sed.End())
        )
      );

    }


    has_deposits_in_tpcset |= (num_deposits_in_tpc > 0);
    std::cout << "Testing TPC " << tpcid << "   " << num_deposits_in_tpc << " (global " << has_deposits_in_tpcset << ")" << std::endl;

  }
  
  std::cout << ">>>> TPCSetID " << tpcset_id << " : " << has_deposits_in_tpcset << std::endl;


  // Get the TPCSet geo of the selected APA/CRM
  // Loop over TPCGeos
  // Get the ActiveBoundingBox geometry
  // Use ContainsPosition() on sim::SimEnergyDeposit start and end point for point with 

  return has_deposits_in_tpcset;
}

DEFINE_ART_MODULE(dunetrigger::HydraFilter)
