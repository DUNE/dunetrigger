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
  //function for matching TPs to sim::IDEs 
  // std::vector<const sim::IDE*> TPToSimIDEs_Ps(recob::Hit const& hit) const;

  //fcl configurable sampling rate 
  // int fADC_SAMPLING_RATE_IN_DTS; // DTS time ticks between adc samples

  // Producer module, configurable from fhicl
  art::InputTag fRawDigitLabel;
  art::InputTag fGenLabel; //generator label for "signal" particles
  std::string fSimChanLabel; // sim channel label

  bool fSaveMC;
  bool fSaveNeutrino;
  bool fSaveRawDigit;

  //ROOT tree for storing output information 
  TTree* fTree; 
  TTree* fRawDigisTree;

  int fEvent, fRun, fSubRun;
  // int fn_TPs; 
  // int fn_noise_TPs_X; 
  // int fn_noise_TPs_U; 
  // int fn_noise_TPs_V; 

  // int fn_signal_TPs_X;
  // int fn_signal_TPs_U;
  // int fn_signal_TPs_V;

  //total charge for the event based on  #electrons at the readout plane
  double total_charge_X;
  double total_charge_U;
  double total_charge_V; 

  //charge of IDEs that were detected via hit finding  for the event 
  // double detected_charge_X;
  // double detected_charge_U;
  // double detected_charge_V;


  //Geant/truth info
  std::vector<int> fNuPDG;
  std::vector<int> fNuCCNC;
  std::vector<int> fNuMode;
  std::vector<float> fNuVx;
  std::vector<float> fNuVy;
  std::vector<float> fNuVz;
  std::vector<float> fNuPx;
  std::vector<float> fNuPy;
  std::vector<float> fNuPz;
  std::vector<float> fNuP;
  std::vector<float> fNuE;


  unsigned int fnParticles; //total number of geant particles 
  std::vector<int>     fTrackId;
  std::vector<float>   fMother; 
  std::vector<float>   fEng;
  std::vector<float>   fEkin;
  std::vector<float>   fMass;
  std::vector<int>     fPdg;
  std::vector<float>   fP;
  std::vector<float>   fPx;
  std::vector<float>   fPy;
  std::vector<float>   fPz;
  std::vector<double>  fstartX;
  std::vector<double>  fstartY;
  std::vector<double>  fstartZ;
  std::vector<double>  fendX;
  std::vector<double>  fendY;
  std::vector<double>  fendZ;


  // raw::RawDigit::ADCvector_t fWaveforms;
  // std::map<raw::ChannelID_t, raw::RawDigit::ADCvector_t> fWaveforms;
  // int fANumber; 


  TH2I* fADCsHistogramByPlane;
  TH1I* fADCsHistogramX;
  TH1I* fADCsHistogramU;
  TH1I* fADCsHistogramV;

  std::map<raw::ChannelID_t, raw::RawDigit::ADCvector_t> fWaveformsBuffer;




};

// -- Constructor --
duneana::RawDigitAna::RawDigitAna(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  fGenLabel(p.get<art::InputTag>("gen_tag")),
  fSimChanLabel(p.get<std::string>("simch_tag", "tpcrawdecoder:simpleSC")),
  fSaveMC(p.get<bool>("SaveMCInfo",true)), 
  fSaveNeutrino(p.get<bool>("SaveNeutrino",false)),
  fSaveRawDigit(p.get<bool>("SaveRawDigits",false))
  {

    if (fSaveRawDigit) {
      fRawDigitLabel = p.get<art::InputTag>("rawdigit_tag"); // Get TP label from fcl file. 

    }

  }


void duneana::RawDigitAna::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs; 
  fTree = tfs->make<TTree>("tree", "Analyser Output Tree");

  fTree->Branch("event",&fEvent,"event/i");
  fTree->Branch("run",&fRun,"run/i");
  fTree->Branch("subrun",&fSubRun,"subrun/i");

  fTree->Branch("totQ_X", &total_charge_X, "totQ_X/d");
  fTree->Branch("totQ_U", &total_charge_U, "totQ_U/d");
  fTree->Branch("totQ_V", &total_charge_V, "totQ_V/d");

  //Nu info
  if (fSaveNeutrino){
    fTree->Branch("NuPDG", &fNuPDG);
    fTree->Branch("NuCCNC", &fNuCCNC);
    fTree->Branch("NuMode", &fNuMode);
    fTree->Branch("NuVx", &fNuVx);
    fTree->Branch("NuVy", &fNuVy);
    fTree->Branch("NuVz", &fNuVz);
    fTree->Branch("NuPx", &fNuPx);
    fTree->Branch("NuPy", &fNuPy);
    fTree->Branch("NuPz", &fNuPz);
    fTree->Branch("NuP", &fNuP);
    fTree->Branch("NuE", &fNuE); 
  }


  //G4 info
  fTree->Branch("nParticles",&fnParticles,"nParticles/i"); 
  fTree->Branch("TrackId",&fTrackId);
  fTree->Branch("Mother",&fMother);
  fTree->Branch("Pdg",&fPdg);
  fTree->Branch("Eng",&fEng);
  fTree->Branch("Ekin",&fEkin);
  fTree->Branch("Mass",&fMass);
  fTree->Branch("P",&fP);
  fTree->Branch("Px",&fPx);
  fTree->Branch("Py",&fPy);
  fTree->Branch("Pz",&fPz);
  fTree->Branch("startX",&fstartX);
  fTree->Branch("startY",&fstartY);
  fTree->Branch("startZ",&fstartZ);
  fTree->Branch("endX",&fendX);
  fTree->Branch("endY",&fendY);
  fTree->Branch("endZ",&fendZ);

  if (fSaveRawDigit){
    std::cout << ">>>>> saving raw digis infos <<<<<" << std::endl;
    fADCsHistogramByPlane = tfs->make<TH2I>("ADCsByPlane", "ADCs by Plane", 4096, 0., 16384., 3, 0., 3.);
    fADCsHistogramX = tfs->make<TH1I>("ADCsPlaneX", "ADCs on plane X", 4096, 0., 16384.);
    fADCsHistogramU = tfs->make<TH1I>("ADCsPlaneU", "ADCs on plane U", 4096, 0., 16384.);
    fADCsHistogramV = tfs->make<TH1I>("ADCsPlaneV", "ADCs on plane X", 4096, 0., 16384.);

    fRawDigisTree = tfs->make<TTree>("rawdigis_tree", "test_ttree");
  
    art::ServiceHandle<geo::Geometry> pgeo;
    auto const& wire_readout = art::ServiceHandle<geo::WireReadout>()->Get();

    for (geo::TPCGeo const& tpc: pgeo->Iterate<geo::TPCGeo>(geo::CryostatID{0})) {

      // Count the channels
      std::vector<raw::ChannelID_t> channels;

      for (auto const& wire : wire_readout.Iterate<geo::WireID>(tpc.ID())) {

          fRawDigisTree->Branch("event",&fEvent,"event/i");
          fRawDigisTree->Branch("run",&fRun,"run/i");
          fRawDigisTree->Branch("subrun",&fSubRun,"subrun/i");

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
  ResetVariables();  // initialise/reset all variables

  fRun = e.run();
  fSubRun = e.subRun();
  fEvent = e.id().event();

  // services 
  art::ServiceHandle<geo::Geometry> geo;
  geo::WireReadoutGeom const &wireReadout = art::ServiceHandle<geo::WireReadout>()->Get();

  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;

  // Load the IDEs for the event to calculate total charge
  auto simchannels = e.getValidHandle<std::vector<sim::SimChannel>>(fSimChanLabel);

  // Process the SimChannels just once
  for (const sim::SimChannel &sc : *simchannels) {
    auto plane =  wireReadout.ROPtoWirePlanes(wireReadout.ChannelToROP(sc.Channel())).at(0).Plane;
    auto tdcidemap = sc.TDCIDEMap();

    // Sum charge for each plane
    for (const auto& tdcide : tdcidemap) {
      for (const auto& ide : tdcide.second) {
        if (plane == geo::kW) {
          total_charge_X += ide.numElectrons;
        }
        else if (plane == geo::kU) {
          total_charge_U += ide.numElectrons;
        }
        else if (plane == geo::kV) {
          total_charge_V += ide.numElectrons;
        }
      }
    }
  }




  if (fSaveRawDigit){

    auto const& digits_handle=e.getValidHandle<std::vector<raw::RawDigit>>(fRawDigitLabel);
    auto& digits_in =*digits_handle;

    for(auto&& digit: digits_in){

      auto plane =  wireReadout.ROPtoWirePlanes(wireReadout.ChannelToROP(digit.Channel())).at(0).Plane;
      int plane_idx = -1;
      TH1I* plane_hist = 0x0;
      switch(plane) {
        case geo::kW:
          plane_idx = 2;
          plane_hist = fADCsHistogramX;
          break;
        case geo::kU:
          plane_hist = fADCsHistogramU;
          plane_idx = 1;
          break;
        case geo::kV:
          plane_hist = fADCsHistogramV;
          plane_idx = 0;
          break;
      }

      for(auto& adc : digit.ADCs() ){
        plane_hist->Fill(adc);
        fADCsHistogramByPlane->Fill(adc,plane_idx);

        fWaveformsBuffer[digit.Channel()]  = digit.ADCs();
      }
    }
  }
  
  //access MC particles through truth label instead of G4 handle, 
  //only interested in "signal" particles &  don't want to iterate over & store info for the avalanche of radiologicals
    
  if (fSaveMC){
    auto mc_truth =  e.getValidHandle< std::vector<simb::MCTruth> >(fGenLabel);
 
    if (mc_truth.isValid()){
      for (const auto& truth : *mc_truth) {

        fnParticles = truth.NParticles();    

        if ( (fSaveNeutrino) && (truth.NeutrinoSet()) ){
          const simb::MCNeutrino& nu = truth.GetNeutrino();
          const simb::MCParticle& neutrino = nu.Nu();

          fNuPDG.push_back(neutrino.PdgCode() ); 
          fNuCCNC.push_back( nu.CCNC() );
          fNuMode.push_back( nu.Mode() );  
          fNuVx.push_back( neutrino.Vx() );
          fNuVy.push_back( neutrino.Vy() );
          fNuVz.push_back( neutrino.Vz() );
          fNuPx.push_back( neutrino.Px() );
          fNuPy.push_back( neutrino.Py() );
          fNuPz.push_back( neutrino.Pz() );
          fNuP.push_back( neutrino.P() );
          fNuE.push_back( neutrino.E() );

        }

        //loop over particles with the signal generator label
        for (int iPartc = 0; iPartc < truth.NParticles(); ++iPartc) {
          
          //grab the current particle
          const simb::MCParticle& trueParticle = truth.GetParticle(iPartc);

          float Ekin = (trueParticle.E() - trueParticle.Mass())*1000;//in MeV
          fTrackId.push_back( trueParticle.TrackId());
          fMother.push_back( trueParticle.Mother());
          fEng.push_back( trueParticle.E() );
          fPdg.push_back( trueParticle.PdgCode());
          fEkin.push_back( Ekin ); 
          fMass.push_back( trueParticle.Mass() );
          fP.push_back( trueParticle.P()) ;
          fPx.push_back( trueParticle.Px());
          fPy.push_back( trueParticle.Py());
          fPz.push_back( trueParticle.Pz());
          fstartX.push_back( trueParticle.Vx());
          fstartY.push_back( trueParticle.Vy());
          fstartZ.push_back( trueParticle.Vz());
          fendX.push_back( trueParticle.EndPosition()[0]);
          fendY.push_back( trueParticle.EndPosition()[1]);
          fendZ.push_back( trueParticle.EndPosition()[2]);
        }
      }
    }
  }
  fTree->Fill();
  if (fSaveRawDigit){
    fRawDigisTree->Fill();
  }

}


void duneana::RawDigitAna::ResetVariables()
{
  fEvent = fRun = fSubRun = -1; 
  // fn_TPs = -1; 

  total_charge_X = total_charge_U = total_charge_V = 0; 

  //Neutrino info 
  fNuPDG.clear();
  fNuCCNC.clear();
  fNuMode.clear();
  fNuVx.clear();
  fNuVy.clear();
  fNuVz.clear();
  fNuPx.clear();
  fNuPy.clear();
  fNuPz.clear();
  fNuP.clear();
  fNuE.clear();


  // GEANT4 info
  fnParticles = 0;
  fTrackId.clear();
  fMother.clear();
  fEng.clear();
  fEkin.clear();
  fMass.clear();
  fPdg.clear();
  fP.clear();
  fPx.clear();
  fPy.clear();
  fPz.clear();
  fstartX.clear();
  fstartY.clear();
  fstartZ.clear();
  fendX.clear();
  fendY.clear();
  fendZ.clear();

  // Waveforms
  fWaveformsBuffer.clear();

}

DEFINE_ART_MODULE(duneana::RawDigitAna)
