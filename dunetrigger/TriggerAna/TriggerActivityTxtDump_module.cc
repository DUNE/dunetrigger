////////////////////////////////////////////////////////////////////////
// Class:       TriggerActivityTxtDump
// Plugin Type: analyzer (Unknown Unknown)
// File:        TriggerActivityTxtDump_module.cc
//
// Generated at Thu Aug 15 14:22:55 2024 by ddrobner using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "detdataformats/trigger/TriggerActivityData.hpp"

#include <fstream>
#include <iostream>

#include <TCanvas.h>
#include <THStack.h>
#include <TLegend.h>
#include <TH1D.h>
#include <art/Framework/Services/Registry/ServiceHandle.h>
#include <fcntl.h>
#include "larcore/Geometry/WireReadout.h"


namespace dunetrigger{
  class TriggerActivityTxtDump;
}

class dunetrigger::TriggerActivityTxtDump : public art::EDAnalyzer {
public:
  explicit TriggerActivityTxtDump(fhicl::ParameterSet const& p);

  // Plugins should not be copied or assigned.
  TriggerActivityTxtDump(TriggerActivityTxtDump const&) = delete;
  TriggerActivityTxtDump(TriggerActivityTxtDump&&) = delete;
  TriggerActivityTxtDump& operator=(TriggerActivityTxtDump const&) = delete;
  TriggerActivityTxtDump& operator=(TriggerActivityTxtDump&&) = delete;

  //void endJob() override;
  void beginJob() override;
  void analyze(art::Event const& e) override;

private:

  art::InputTag sim_ta_tag;
  std::ofstream sim_out;
  geo::WireReadoutGeom const *geom_;


};


dunetrigger::TriggerActivityTxtDump::TriggerActivityTxtDump(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  sim_ta_tag(p.get<art::InputTag>("sim_ta_tag")),
  geom_(&art::ServiceHandle<geo::WireReadout>()->Get())
    // ,
  // More initializers here.
{
  using dunedaq::trgdataformats::TriggerActivityData;
  // Call appropriate consumes<>() for any products to be retrieved by this
  // module.
  consumes<std::vector<TriggerActivityData>>(sim_ta_tag);
}

void dunetrigger::TriggerActivityTxtDump::beginJob() {
  sim_out.open("sim_tas.txt");
}

void dunetrigger::TriggerActivityTxtDump::analyze(art::Event const& e)
{
  using dunedaq::trgdataformats::TriggerActivityData;
  // Implementation of required member function here.
  std::vector<TriggerActivityData> sim_tas = *(e.getValidHandle<std::vector<TriggerActivityData>>(sim_ta_tag));

  std::cout << "Sim Length: " << sim_tas.size() << std::endl;

  

  // order is time_start, time_peak, time_activity, adc_peak, channel_peak,
  // channel_start, adc_int, event, time_end
  for (auto s : sim_tas) {
    readout::ROPID rop    = geom_->ChannelToROP(s.channel_start);
    readout::TPCsetID tpcset = rop.asTPCsetID();

    sim_out << e.event()       << ","
            << e.run()         << ","
            << e.subRun()      << ","
            << s.time_start    << ","
            << s.time_peak     << ","
            << s.time_end      << ","
            << s.time_activity << ","
            << s.channel_start << ","
            << s.channel_peak  << ","
            << s.channel_end   << ","
            << s.adc_integral  << ","
            << s.adc_peak      << ","
            << (uint16_t)s.algorithm << ","
            << (uint16_t)s.detid     << ","
            << (uint16_t)s.type      << ","
            << rop.ROP               << ","
            << tpcset.TPCset         << std::endl;
  }

}

DEFINE_ART_MODULE(dunetrigger::TriggerActivityTxtDump)


