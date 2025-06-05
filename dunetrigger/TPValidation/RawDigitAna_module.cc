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
  std::vector<const sim::IDE*> TPToSimIDEs_Ps(recob::Hit const& hit) const;

  //fcl configurable sampling rate 
  int fADC_SAMPLING_RATE_IN_DTS; // DTS time ticks between adc samples

  // Producer module, configurable from fhicl
  art::InputTag fRawDigitLabel;
  art::InputTag fTPLabel;
  art::InputTag fGenLabel; //generator label for "signal" particles
  std::string fSimChanLabel; // sim channel label

  bool fSaveMC;
  bool fSaveNeutrino;
  bool fSaveRawDigit;
  bool fSaveTPs; 

  int foffsetU; //offset for valid TP window for IDE matching (needed for induction sigs).  
  int foffsetV; 

  //ROOT tree for storing output information 
  TTree* fTree; 
  // TTree* fRawDigisTree;

  int fEvent, fRun, fSubRun;
  int fn_TPs; 
  int fn_noise_TPs_X; 
  int fn_noise_TPs_U; 
  int fn_noise_TPs_V; 

  int fn_signal_TPs_X;
  int fn_signal_TPs_U;
  int fn_signal_TPs_V;

  //total charge for the event based on  #electrons at the readout plane
  double total_charge_X;
  double total_charge_U;
  double total_charge_V; 

  //charge of IDEs that were detected via hit finding  for the event 
  double detected_charge_X;
  double detected_charge_U;
  double detected_charge_V;


  // Vectors for storing Trigger Primitive (TP) data 
  std::vector<uint32_t> fTP_channels;
  std::vector<uint64_t> fTP_startT; 
  std::vector<uint64_t> fTP_peakT; 
  std::vector<uint64_t> fTP_TOT; 
  std::vector<uint32_t> fTP_SADC;
  std::vector<uint16_t> fTP_peakADC; 
  std::vector<int>      fTP_plane;
  std::vector<int>      fTP_TPC;
  std::vector<int>      fTP_trueX;
  std::vector<int>      fTP_trueY;
  std::vector<int>      fTP_trueZ;
  std::vector<int>      fTP_signal;


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




};

// -- Constructor --
duneana::RawDigitAna::RawDigitAna(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  fADC_SAMPLING_RATE_IN_DTS(p.get<int>("ADC_SAMPLING_RATE_IN_DTS",32)),
  fGenLabel(p.get<art::InputTag>("gen_tag")),
  fSimChanLabel(p.get<std::string>("simch_tag", "tpcrawdecoder:simpleSC")),
  // fRawDigitLabel(p.get<std::string>("rawdigit_tag", "tpcrawdecoder:daq")),
  fSaveMC(p.get<bool>("SaveMCInfo",true)), 
  fSaveNeutrino(p.get<bool>("SaveNeutrino",false)),
  fSaveRawDigit(p.get<bool>("SaveRawDigits",false)),
  fSaveTPs(p.get<bool>("SaveTPInfo",true)),
  foffsetU(p.get<int>("U_window_end_offset",11)),
  foffsetV(p.get<int>("V_window_end_offset",10))
  {

    if (fSaveTPs) {
      fTPLabel = p.get<art::InputTag>("tp_tag"); // Get TP label from fcl file. 
    }
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

  // TPG info
  fTree->Branch("n_TPs", &fn_TPs, "n_TPs/i");
  fTree->Branch("n_noise_TPs_X", &fn_noise_TPs_X, "n_noise_TPs_X/i");
  fTree->Branch("n_noise_TPs_U", &fn_noise_TPs_U, "n_noise_TPs_U/i");
  fTree->Branch("n_noise_TPs_V", &fn_noise_TPs_V, "n_noise_TPs_V/i");

  fTree->Branch("n_signal_TPs_X", &fn_signal_TPs_X, "n_signal_TPs_X/i");
  fTree->Branch("n_signal_TPs_U", &fn_signal_TPs_U, "n_signal_TPs_U/i");
  fTree->Branch("n_signal_TPs_V", &fn_signal_TPs_V, "n_signal_TPs_V/i");  


  fTree->Branch("totQ_X", &total_charge_X, "totQ_X/d");
  fTree->Branch("totQ_U", &total_charge_U, "totQ_U/d");
  fTree->Branch("totQ_V", &total_charge_V, "totQ_V/d");

  fTree->Branch("detQ_X", &detected_charge_X, "detQ_X/d");
  fTree->Branch("detQ_U", &detected_charge_U, "detQ_U/d");
  fTree->Branch("detQ_V", &detected_charge_V, "detQ_V/d");


  //TP info
  if (fSaveTPs){ 
    fTree->Branch("TP_channel", &fTP_channels);
    fTree->Branch("TP_startT", &fTP_startT);
    fTree->Branch("TP_peakT", &fTP_peakT);
    fTree->Branch("TP_TOT", &fTP_TOT);
    fTree->Branch("TP_SADC", &fTP_SADC);
    fTree->Branch("TP_peakADC", &fTP_peakADC);
    fTree->Branch("TP_plane", &fTP_plane);
    fTree->Branch("TP_TPC", &fTP_TPC);
    fTree->Branch("TP_trueX", &fTP_trueX);
    fTree->Branch("TP_trueY", &fTP_trueY);
    fTree->Branch("TP_trueZ", &fTP_trueZ);
    fTree->Branch("TP_signal", &fTP_signal);  
  }


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
    // fRawDigisTree = tfs->make<TTree>("raw_digis", "Analyser Output Tree");
    // fADCsHistogramByPlane = new TH2I("ADCsByPlane", "ADCs by Plane", 4096, 0., 16384., 3, 0., 2.);
    // fRawDigisTree->Branch("adcs_by_plane", &fADCsHistogramByPlane); 
    // fRawDigisTree->Branch("a_number", &fANumber); 
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
  //auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e); 

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
          plane_idx = 1;
          plane_hist = fADCsHistogramX;
          break;
        case geo::kU:
          plane_hist = fADCsHistogramU;
          plane_idx = 2;
          break;
        case geo::kV:
          plane_hist = fADCsHistogramV;
          plane_idx = 3;
          break;
      }

      for(auto& adc : digit.ADCs() ){
        plane_hist->Fill(adc);
        fADCsHistogramByPlane->Fill(adc,plane_idx);
      }
      // fWaveforms.insert(fWaveforms.end(), digit.ADCs().begin(), digit.ADCs().end());
      
      
      // std::cout << "Channel: " << digit.Channel() << " plane " << plane << " plane-idx " << plane_idx << "  first ADC: " <<  digit.ADCs()[0] << std::endl;
      //  << std::endl;
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
  // if (fSaveRawDigit){
    // fRawDigisTree->Fill();
  // }

}


//TP - ide matching based on hit start and end times rather than peak ADC time 
std::vector<const sim::IDE*> duneana::RawDigitAna::TPToSimIDEs_Ps(recob::Hit const& hit) const
{
  std::vector<const sim::IDE*> retVec;
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  
  //IDEs are computed based on the time during which the ionization cloud hits the readout plane. For ind. signals this corresponds to the inflection point
  //need to introduce an offset for this effect to be properly accounted for 
  int offset = 0;

  if (hit.View() == geo::kU) offset = foffsetU; //fcl configurable offset for induction TPs 
  else if (hit.View() == geo::kV) offset = foffsetV;  

  // Adjust the start and end time based on the hit, this may need to be expanded or modified
  int start_tdc = hit.StartTick(); //clockData.TPCTick2TDC(hit.PeakTimeMinusRMS(fHitTimeRMS));
  int end_tdc = hit.EndTick() + offset;//clockData.TPCTick2TDC(hit.PeakTimePlusRMS(fHitTimeRMS));

  // Ensure no negative TDC values
  if (start_tdc < 0) start_tdc = 0;
  if (end_tdc < 0) end_tdc = 0;

  if (start_tdc > end_tdc) { throw; }

  // Get the map of TDC values to IDEs
  const std::vector<std::pair<unsigned short, std::vector<sim::IDE>>>& tdcIDEMap = (bt_serv->FindSimChannel(hit.Channel()))->TDCIDEMap();

  // Create a vector of pointers to the elements in the TDC-IDE map
  std::vector<const std::pair<unsigned short, std::vector<sim::IDE>>*> tdcIDEMap_SortedPointers;
  for (auto& pair : tdcIDEMap) {
    tdcIDEMap_SortedPointers.push_back(&pair);
  }

  // Sort the map based on TDC values
  auto pairSort = [](auto& a, auto& b) { return a->first < b->first; };
  if (!std::is_sorted(tdcIDEMap_SortedPointers.begin(), tdcIDEMap_SortedPointers.end(), pairSort)) {
    std::sort(tdcIDEMap_SortedPointers.begin(), tdcIDEMap_SortedPointers.end(), pairSort);
  }

  // Create dummy vectors for comparing pairs
  std::vector<sim::IDE> dummyVec;
  std::pair<double, std::vector<sim::IDE>> start_tdcPair = std::make_pair(start_tdc, dummyVec);
  std::pair<double, std::vector<sim::IDE>> end_tdcPair = std::make_pair(end_tdc, dummyVec);

  auto start_tdcPair_P = &start_tdcPair;
  auto end_tdcPair_P = &end_tdcPair;

  // Find the range of IDEs that fall between start_tdc and end_tdc using lower and upper bounds
  auto mapFirst = std::lower_bound(tdcIDEMap_SortedPointers.begin(), tdcIDEMap_SortedPointers.end(), start_tdcPair_P, pairSort);

  // Upper bound is exclusive --> should give the first element after the end_tdc
  auto mapLast = std::upper_bound(tdcIDEMap_SortedPointers.begin(), tdcIDEMap_SortedPointers.end(), end_tdcPair_P, pairSort);

  // Iterate through the range of IDEs that fall between start_tdc and end_tdc and save them to a vector
  for (auto& mapitr = mapFirst; mapitr != mapLast; ++mapitr) {
    for (auto& ide : (*mapitr)->second) {
      // Save all IDEs within the specified time window
      retVec.push_back(&ide);
    }
  }
  // output the IDE vector 
  return retVec;
}



void duneana::RawDigitAna::ResetVariables()
{
  fEvent = fRun = fSubRun = -1; 
  fn_TPs = -1; 


  total_charge_X = total_charge_U = total_charge_V = 0; 


  detected_charge_X = detected_charge_U = detected_charge_V = 0; 

  //initialise number of noise and signal TPs for this event 
  fn_noise_TPs_X = fn_noise_TPs_U = fn_noise_TPs_V =  0;
  fn_signal_TPs_X = fn_signal_TPs_U = fn_signal_TPs_V =  0;

  // fADCsHistogramByPlane->Clear();

  fTP_channels.clear(); 
  fTP_startT.clear();
  fTP_peakT.clear();
  fTP_TOT.clear();
  fTP_SADC.clear();
  fTP_peakADC.clear();
  fTP_plane.clear();
  fTP_TPC.clear();
  fTP_trueX.clear();
  fTP_trueY.clear();
  fTP_trueZ.clear();
  fTP_signal.clear();


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

}

DEFINE_ART_MODULE(duneana::RawDigitAna)
