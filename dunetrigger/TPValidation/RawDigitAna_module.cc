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

  void reset_variables();

  // Producer module, configurable from fhicl
  art::InputTag fRawDigitLabel; // raw digi label
  std::string fSimChanLabel; // sim channel label

  bool fSaveWaveforms;
  bool fSaveChannelsWithIDEsOnly;

  //ROOT tree for storing output information 
  TTree* fRawDigisTree = nullptr;

  int fEvent, fRun, fSubRun;

  TH2I* fADCsHistogramByPlane = nullptr;
  TH1I* fADCsHistogramX = nullptr;
  TH1I* fADCsHistogramU = nullptr;
  TH1I* fADCsHistogramV = nullptr;

  TH2I* fNoiseADCsHistogramByPlane = nullptr;
  TH1I* fNoiseADCsHistogramX = nullptr;
  TH1I* fNoiseADCsHistogramU = nullptr;
  TH1I* fNoiseADCsHistogramV = nullptr;

  std::map<raw::ChannelID_t, raw::RawDigit::ADCvector_t> fWaveformsBuffer;
  std::vector<unsigned int> fActiveChannels;

};

// -- Constructor --
duneana::RawDigitAna::RawDigitAna(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  fRawDigitLabel(p.get<art::InputTag>("rawdigit_tag", "tpcrawdecoder:daq")),
  fSimChanLabel(p.get<std::string>("simch_tag", "tpcrawdecoder:simpleSC")),
  fSaveWaveforms(p.get<bool>("save_waveforms",true)),
  fSaveChannelsWithIDEsOnly(p.get<bool>("save_channels_with_ides_only",false))
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

  fNoiseADCsHistogramByPlane = tfs->make<TH2I>("ADCsNoiseByPlane", "ADCs by Plane", 4096, 0., 16384., 3, 0., 3.);
  fNoiseADCsHistogramX = tfs->make<TH1I>("ADCsNoisePlaneX", "ADCs on plane X", 4096, 0., 16384.);
  fNoiseADCsHistogramU = tfs->make<TH1I>("ADCsNoisePlaneU", "ADCs on plane U", 4096, 0., 16384.);
  fNoiseADCsHistogramV = tfs->make<TH1I>("ADCsNoisePlaneV", "ADCs on plane X", 4096, 0., 16384.);


  fRawDigisTree = tfs->make<TTree>("rawdigis_tree", "test_ttree");

  fRawDigisTree->Branch("event",&fEvent,"event/i");
  fRawDigisTree->Branch("run",&fRun,"run/i");
  fRawDigisTree->Branch("subrun",&fSubRun,"subrun/i");

  if (fSaveChannelsWithIDEsOnly) {
    fRawDigisTree->Branch("active_channels", &fActiveChannels);
  }

  if ( fSaveWaveforms ) {

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


}
void duneana::RawDigitAna::analyze(art::Event const& e)
{
  reset_variables();  // initialise/reset all variables

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

  fActiveChannels.reserve(signal_chans.size());  // optional but avoids reallocations
  std::copy(signal_chans.begin(), signal_chans.end(), std::back_inserter(fActiveChannels));

  auto const& digits_handle=e.getValidHandle<std::vector<raw::RawDigit>>(fRawDigitLabel);
  auto& digits_in =*digits_handle;

  size_t channel_count(0);
  size_t digi_count(0);

  std::cout << digits_in.size() << " Digits found in this event" << std::endl;

  for(auto&& digit: digits_in){

    if ( channel_count % 1000 == 0) 
      std::cout << " >>> " << channel_count << std::endl;
    ++channel_count;

    auto plane =  wireReadout.ROPtoWirePlanes(wireReadout.ChannelToROP(digit.Channel())).at(0).Plane;
    int plane_idx = -1;
    TH1I* plane_adc_hist(0x0);
    TH1I* plane_noise_hist(0x0);
    switch(plane) {

      case geo::kU:
        plane_adc_hist = fADCsHistogramU;
        plane_noise_hist = fNoiseADCsHistogramU;
        plane_idx = 0;
        break;
      case geo::kV:
        plane_adc_hist = fADCsHistogramV;
        plane_noise_hist = fNoiseADCsHistogramV;
        plane_idx = 1;
        break;
      case geo::kW:
        plane_adc_hist = fADCsHistogramX;
        plane_noise_hist = fNoiseADCsHistogramX;
        plane_idx = 2;
        break;
    }

    bool has_signal = (signal_chans.find(digit.Channel()) != signal_chans.end());
    // Fill histograms
    for(auto& adc : digit.ADCs() ){

      plane_adc_hist->Fill(adc);
      fADCsHistogramByPlane->Fill(adc,plane_idx);

      // Fill noise histogram with adcs where there are no IDEs
      if (!has_signal) {
        plane_noise_hist->Fill(adc);
        fNoiseADCsHistogramByPlane->Fill(adc,plane_idx);
      }

    }

    // Save ADCs if they have signals, or only save signals wf is disabled
    if ( fSaveWaveforms and (has_signal or !fSaveChannelsWithIDEsOnly) ) {
      fWaveformsBuffer[digit.Channel()]  = digit.ADCs();        
    }

    digi_count += digit.ADCs().size();

  }

  std::cout << "Processed channels " << channel_count << std::endl;
  std::cout << "Processed rawdigit samples " << digi_count << std::endl;

  if (fSaveWaveforms)
    fRawDigisTree->Fill();

}


void duneana::RawDigitAna::reset_variables()
{
  // Run information
  fEvent = fRun = fSubRun = -1; 
  fActiveChannels.clear();
  // Waveforms
  for( auto& [ch, buffer] : fWaveformsBuffer  ) {
    buffer.clear();
  }

}

DEFINE_ART_MODULE(duneana::RawDigitAna)
