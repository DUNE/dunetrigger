////////////////////////////////////////////////////////////////////////
// Class:       RawDigitsCleaner
// Plugin Type: producer
// File:        RawDigitsCleaner_module.cc
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/WireReadout.h"

#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/Simulation/SimChannel.h"

#include <memory>
#include <regex>

namespace dunetrigger {
class RawDigitsCleaner;
}

class dunetrigger::RawDigitsCleaner : public art::EDProducer {
public:
  explicit RawDigitsCleaner(fhicl::ParameterSet const& p);

  RawDigitsCleaner(RawDigitsCleaner const&) = delete;
  RawDigitsCleaner(RawDigitsCleaner&&) = delete;
  RawDigitsCleaner& operator=(RawDigitsCleaner const&) = delete;
  RawDigitsCleaner& operator=(RawDigitsCleaner&&) = delete;

  void produce(art::Event& e) override;

private:

  // TODO: rename n_tpsets
  int num_tpcset;

  std::string get_prodid(std::string prod_base, int i) {
    // C++20 where are you?
    std::ostringstream oss;
    oss << prod_base << i;
    return oss.str();
    // return prod_base+std::to_string(i);
  }
};

dunetrigger::RawDigitsCleaner::RawDigitsCleaner(fhicl::ParameterSet const& p)
    : EDProducer{p}
{

  // get a service handle for geometry
  geo::WireReadoutGeom const *wrgeom = &art::ServiceHandle<geo::WireReadout>()->Get();
  
  int num_tpcset = wrgeom->MaxTPCsets();
  for(int i(0); i<num_tpcset; ++i) {

    // Move to config parameters
    std::string rawdigit_name = get_prodid("daq", i);
    std::string simchan_name = get_prodid("simpleSC", i);

    // std::cout << "Registering '" << rawdigit_name << "'  '" << simchan_name << "'" << std::endl;
    produces<std::vector<raw::RawDigit>>(rawdigit_name);
    produces<std::vector<sim::SimChannel>>(simchan_name);
  }
}

void dunetrigger::RawDigitsCleaner::produce(art::Event& e)
{


    // TODO: alternative implementation that does not rely on a config parameter
    // Get the number of TPCSets from the wiregeometry
    // Loop from 0 to NTPCSets
    // getValidHandle("simpleSC{i_tpcset}")
    // if doesn't exist -> handle
    // else continue as it is
  for(int i(0); i<num_tpcset; ++i) {

    std::string rawdigit_name = get_prodid("daq", i);
    std::string simchan_name = get_prodid("simpleSC", i);

    // TODO: check for handle validity
    auto rawdigit_handle = e.getValidHandle<std::vector<raw::RawDigit>>("tpcrawdecoder:"+rawdigit_name);
    auto simchan_handle = e.getValidHandle<std::vector<sim::SimChannel>>("tpcrawdecoder:"+simchan_name);

    size_t n_chans_with_activity{0};

    for( auto sc : *simchan_handle) {
      if (!sc.TDCIDEMap().empty()) ++n_chans_with_activity;
    }

    std::cout << simchan_name << "  " << n_chans_with_activity << std::endl;

    auto digits = std::make_unique<std::vector<raw::RawDigit>>();
    auto simchans = std::make_unique<std::vector<sim::SimChannel>>();
    if (n_chans_with_activity > 0) {
      *digits = *rawdigit_handle;
      *simchans = *simchan_handle;
    }
    e.put(std::move(digits), rawdigit_name);
    e.put(std::move(simchans), simchan_name);
  }
}

DEFINE_ART_MODULE(dunetrigger::RawDigitsCleaner)
