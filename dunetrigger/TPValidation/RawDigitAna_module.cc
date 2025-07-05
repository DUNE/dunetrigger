////////////////////////////////////////////////////////////////////////
// Class:       RawDigitAna
// Plugin Type: analyzer
// File:        RawDigitAna_module.cc
//
// Generated at Fri Nov 15 01:56:30 2024 by Klaudia Wawrowska
//
// Module for dumping TPG and associated MC information to TTree for analysis.
//
// Simulation producers to specify: tp_tag, gen_tag (generator - singles, marley, genie etc.), simch_tag, 
// Booleans which allow user to choose which information to save: SaveNeutrino, SaveMC (MCTruth not G4), SaveTPs
// 
// The association between MC and TP information uses IDEs. This requires some tweaks to make it work 
// Specify ADC_SAMPLING_RATE_IN_DTS to go back to TPC ticks used in the simulation (1 tick = 500 ns or whatever) 
// U/V window offsets for algs which work only on +ve induction pulses has to be set for 11 for U and 10/9 for V.
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


#include "detdataformats/trigger/TriggerCandidateData.hpp"
#include "detdataformats/trigger/TriggerActivityData.hpp"
#include "detdataformats/trigger/TriggerPrimitive.hpp"
#include "detdataformats/DetID.hpp"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/RDTimeStamp.h"
//#include "lardata/DetectorInfoServices/DetectorClocksService.h"


#include "larcore/Geometry/Geometry.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"

// Simulation truth
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/MCBase/MCStep.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larcore/Geometry/WireReadout.h"
#include "lardataobj/RecoBase/Hit.h"

// ROOT includes
#include <TH1I.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH2I.h>
#include <TTree.h>
#include <TFile.h>
#include <TGraph.h>
#include <TRandom.h>
#include <TVector3.h>
#include <fcntl.h>

#include <memory>
#include <algorithm>
#include <iostream>



namespace duneana {
  class RawDigitAna;
}


class duneana::RawDigitAna : public art::EDAnalyzer {
public:
  explicit RawDigitAna(fhicl::ParameterSet const& p);

  RawDigitAna(RawDigitAna const&) = delete;
  RawDigitAna(RawDigitAna&&) = delete;
  RawDigitAna& operator=(RawDigitAna const&) = delete;
  RawDigitAna& operator=(RawDigitAna&&) = delete;

  
  void analyze(art::Event const& e) override;
  void beginJob() override;
 
private:

  void ResetVariables();

  // Producer module, configurable from fhicl
  art::InputTag fRawDigitLabel; // raw digi label
  std::string fSimChanLabel; // sim channel label

  bool fSaveSignalWaveformsOnly;

  //ROOT tree for storing output information 
  TTree* fRawDigisTree;

  int fEvent, fRun, fSubRun;

  TH2I* fADCsHistogramByPlane;
  TH1I* fADCsHistogramX;
  TH1I* fADCsHistogramU;
  TH1I* fADCsHistogramV;

  std::map<raw::ChannelID_t, raw::RawDigit::ADCvector_t> fWaveformsBuffer;
  std::vector<unsigned int> fChannelsWithelectrons;

};

// -- Constructor --
duneana::RawDigitAna::RawDigitAna(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  fRawDigitLabel(p.get<art::InputTag>("rawdigit_tag", "tpcrawdecoder:daq")),
  fSimChanLabel(p.get<std::string>("simch_tag", "tpcrawdecoder:simpleSC")),
  fSaveSignalWaveformsOnly(p.get<bool>("SaveSignalOnly",false))
  {
  }


void duneana::RawDigitAna::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs; 

  std::cout << ">>>>> saving raw digis infos <<<<<" << std::endl;
  fADCsHistogramByPlane = tfs->make<TH2I>("ADCsByPlane", "ADCs by Plane", 4096, 0., 16384., 3, 0., 3.);
  fADCsHistogramX = tfs->make<TH1I>("ADCsPlaneX", "ADCs on plane X", 4096, 0., 16384.);
  fADCsHistogramU = tfs->make<TH1I>("ADCsPlaneU", "ADCs on plane U", 4096, 0., 16384.);
  fADCsHistogramV = tfs->make<TH1I>("ADCsPlaneV", "ADCs on plane X", 4096, 0., 16384.);

  fRawDigisTree = tfs->make<TTree>("rawdigis_tree", "test_ttree");

  fRawDigisTree->Branch("event",&fEvent,"event/i");
  fRawDigisTree->Branch("run",&fRun,"run/i");
  fRawDigisTree->Branch("subrun",&fSubRun,"subrun/i");

  if (fSaveSignalWaveformsOnly) {
    fRawDigisTree->Branch("chans_with_electrons", &fChannelsWithelectrons);
  }

  art::ServiceHandle<geo::Geometry> pgeo;
  auto const& wire_readout = art::ServiceHandle<geo::WireReadout>()->Get();

  for (geo::TPCGeo const& tpc: pgeo->Iterate<geo::TPCGeo>(geo::CryostatID{0})) {

    // Count the channels
    std::vector<raw::ChannelID_t> channels;

    for (auto const& wire : wire_readout.Iterate<geo::WireID>(tpc.ID())) {

        raw::ChannelID_t ch = wire_readout.PlaneWireToChannel(wire);
        channels.push_back(ch);
        
        auto& buffer = fWaveformsBuffer[ch];
        fRawDigisTree->Branch(std::to_string(ch).c_str(), &buffer);
    }
    std::cout << "N Channels: " << channels.size() << "(First channel " << channels[0] << ")" << std::endl;
  }


}
void duneana::RawDigitAna::analyze(art::Event const& e)
{
  ResetVariables();  // initialise/reset all variables

  fRun = e.run();
  fSubRun = e.subRun();
  fEvent = e.id().event();

  // services 
  art::ServiceHandle<geo::Geometry> geo;
  geo::WireReadoutGeom const &wireReadout = art::ServiceHandle<geo::WireReadout>()->Get();

  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;

  auto simchans_handle = e.getValidHandle<std::vector<sim::SimChannel>>(fSimChanLabel);
  
  std::set<unsigned int> signal_chans;
  for (const sim::SimChannel &sc : *simchans_handle) {

    if (!sc.TDCIDEMap().empty())
      signal_chans.insert(sc.Channel());
  }
  std::cout << signal_chans.size() << " channels with ides found in this event" << std::endl;

  fChannelsWithelectrons.reserve(signal_chans.size());  // optional but avoids reallocations
  std::copy(signal_chans.begin(), signal_chans.end(), std::back_inserter(fChannelsWithelectrons));

  auto const& digits_handle=e.getValidHandle<std::vector<raw::RawDigit>>(fRawDigitLabel);
  auto& digits_in =*digits_handle;

  size_t counts(0);

  std::cout << digits_in.size() << " Digits found in this event" << std::endl;

  for(auto&& digit: digits_in){

    if ( counts % 1000 == 0) 
      std::cout << " >>> " << counts << std::endl;
    ++counts;

    if ( fSaveSignalWaveformsOnly and !(signal_chans.find(digit.Channel()) != signal_chans.end()))
      continue;

    auto plane =  wireReadout.ROPtoWirePlanes(wireReadout.ChannelToROP(digit.Channel())).at(0).Plane;
    int plane_idx = -1;
    TH1I* plane_hist = 0x0;
    switch(plane) {

      case geo::kU:
        plane_hist = fADCsHistogramU;
        plane_idx = 0;
        break;
      case geo::kV:
        plane_hist = fADCsHistogramV;
        plane_idx = 1;
        break;
      case geo::kW:
        plane_idx = 2;
        plane_hist = fADCsHistogramX;
        break;
    }

    for(auto& adc : digit.ADCs() ){
      plane_hist->Fill(adc);
      fADCsHistogramByPlane->Fill(adc,plane_idx);

      fWaveformsBuffer[digit.Channel()]  = digit.ADCs();
    }
  }

  fRawDigisTree->Fill();

}


void duneana::RawDigitAna::ResetVariables()
{
  // Run information
  fEvent = fRun = fSubRun = -1; 
  fChannelsWithelectrons.clear();
  // Waveforms
  for( auto& [ch, buffer] : fWaveformsBuffer  ) {
    buffer.clear();
  }

}

DEFINE_ART_MODULE(duneana::RawDigitAna)
