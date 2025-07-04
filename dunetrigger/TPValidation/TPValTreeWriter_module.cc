////////////////////////////////////////////////////////////////////////
// Class:       TPValTreeWriter
// Plugin Type: analyzer
// File:        TPGAna_module.cc
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
#include <TNamed.h>
#include <TVector3.h>
#include <fcntl.h>

#include <memory>
#include <algorithm>
#include <iostream>

#include <nlohmann/json.hpp>
using json = nlohmann::json;


namespace duneana {


struct EventDataBuffer {
  int event;
  int run;
  int subrun;

  void branch_on( TTree* tree ) {
    tree->Branch("event",&event,"event/i");
    tree->Branch("run",&run,"run/i");
    tree->Branch("subrun",&subrun,"subrun/i");
  }

  void clear() {
    event = -1;
    run = -1;
    subrun = -1; 
  }
};

struct MCTruthDataBuffer {
  // Monte carlo truth information
  unsigned int        num_particles; //total number of geant particles 
  std::vector<int>    TrackId;
  std::vector<float>  Mother; 
  std::vector<float>  Eng;
  std::vector<float>  Ekin;
  std::vector<float>  Mass;
  std::vector<int>    Pdg;
  std::vector<float>  P;
  std::vector<float>  Px;
  std::vector<float>  Py;
  std::vector<float>  Pz;
  std::vector<double> startX;
  std::vector<double> startY;
  std::vector<double> startZ;
  std::vector<double> endX;
  std::vector<double> endY;
  std::vector<double> endZ;

  void branch_on( TTree* tree ) {
    tree->Branch("nParticles",&num_particles,"nParticles/i"); 
    tree->Branch("TrackId",&TrackId);
    tree->Branch("Mother",&Mother);
    tree->Branch("Pdg",&Pdg);
    tree->Branch("Eng",&Eng);
    tree->Branch("Ekin",&Ekin);
    tree->Branch("Mass",&Mass);
    tree->Branch("P",&P);
    tree->Branch("Px",&Px);
    tree->Branch("Py",&Py);
    tree->Branch("Pz",&Pz);
    tree->Branch("startX",&startX);
    tree->Branch("startY",&startY);
    tree->Branch("startZ",&startZ);
    tree->Branch("endX",&endX);
    tree->Branch("endY",&endY);
    tree->Branch("endZ",&endZ); 
  }

  void clear() {
   num_particles = 0;
   TrackId.clear();
   Mother.clear();
   Eng.clear();
   Ekin.clear();
   Mass.clear();
   Pdg.clear();
   P.clear();
   Px.clear();
   Py.clear();
   Pz.clear();
   startX.clear();
   startY.clear();
   startZ.clear();
   endX.clear();
   endY.clear();
   endZ.clear();
  }

};

struct IDEsDataBuffer {

  // Data buffers
  unsigned int entries; //total number of geant particles 
  std::vector<uint32_t> channel;
  std::vector<float> time;
  std::vector<int> track_id;
  std::vector<float> numElectrons;
  std::vector<float> energy;
  std::vector<float> x;
  std::vector<float> y;
  std::vector<float> z;

  void branch_on( TTree* tree ) {
    tree->Branch("nIDEs",&entries,"nIDEs/i");
    tree->Branch("channel",&channel);  
    tree->Branch("time",&time);  
    tree->Branch("TrackId",&track_id);  
    tree->Branch("nElectrons",&numElectrons);  
    tree->Branch("energy",&energy);  
    tree->Branch("x",&x);  
    tree->Branch("y",&y);  
    tree->Branch("z",&z);  
  }
  
  void clear() {
    entries = 0;
    channel.clear();
    time.clear();
    track_id.clear();
    numElectrons.clear();
    energy.clear();
    x.clear();
    y.clear();
    z.clear();
  }
};

struct TPSummaryDataBuffer {

  int num_TPs; 
  int num_noise_TPs_X; 
  int num_noise_TPs_U; 
  int num_noise_TPs_V; 

  int num_signal_TPs_X;
  int num_signal_TPs_U;
  int num_signal_TPs_V;

  //total charge for the event based on  #electrons at the readout plane
  double total_charge_X;
  double total_charge_U;
  double total_charge_V; 

  //charge of IDEs that were detected via hit finding  for the event 
  double detected_charge_X;
  double detected_charge_U;
  double detected_charge_V;

  void branch_on( TTree* tree ) {
    // TPG info
    tree->Branch("n_TPs", &num_TPs, "n_TPs/i");
    tree->Branch("n_noise_TPs_U", &num_noise_TPs_U, "n_noise_TPs_U/i");
    tree->Branch("n_noise_TPs_V", &num_noise_TPs_V, "n_noise_TPs_V/i");
    tree->Branch("n_noise_TPs_X", &num_noise_TPs_X, "n_noise_TPs_X/i");

    tree->Branch("n_signal_TPs_U", &num_signal_TPs_U, "n_signal_TPs_U/i");
    tree->Branch("n_signal_TPs_V", &num_signal_TPs_V, "n_signal_TPs_V/i");  
    tree->Branch("n_signal_TPs_X", &num_signal_TPs_X, "n_signal_TPs_X/i");


    tree->Branch("totQ_U", &total_charge_U, "totQ_U/d");
    tree->Branch("totQ_V", &total_charge_V, "totQ_V/d");
    tree->Branch("totQ_X", &total_charge_X, "totQ_X/d");
    
    tree->Branch("detQ_U", &detected_charge_U, "detQ_U/d");
    tree->Branch("detQ_V", &detected_charge_V, "detQ_V/d");
    tree->Branch("detQ_X", &detected_charge_X, "detQ_X/d");
  }

  void clear() {
    num_TPs = -1; 

    total_charge_X = total_charge_U = total_charge_V = 0; 

    detected_charge_X = detected_charge_U = detected_charge_V = 0; 

    //initialise number of noise and signal TPs for this event 
    num_noise_TPs_X = num_noise_TPs_U = num_noise_TPs_V =  0;
    num_signal_TPs_X = num_signal_TPs_U = num_signal_TPs_V =  0;
  }

};

struct TPsDataBuffer {

  // Vectors for storing Trigger Primitive (TP) data 
  std::vector<uint32_t> channels;
  std::vector<uint64_t> startT; 
  std::vector<uint64_t> peakT; 
  std::vector<uint64_t> TOT; 
  std::vector<uint32_t> SADC;
  std::vector<uint16_t> peakADC; 
  std::vector<int>      plane;
  std::vector<int>      TPC;
  std::vector<int>      trueX;
  std::vector<int>      trueY;
  std::vector<int>      trueZ;
  std::vector<int>      signal;

  void branch_on( TTree* tree ) {
    tree->Branch("TP_channel", &channels);
    tree->Branch("TP_startT", &startT);
    tree->Branch("TP_peakT", &peakT);
    tree->Branch("TP_TOT", &TOT);
    tree->Branch("TP_SADC", &SADC);
    tree->Branch("TP_peakADC", &peakADC);
    tree->Branch("TP_plane", &plane);
    tree->Branch("TP_TPC", &TPC);
    tree->Branch("TP_trueX", &trueX);
    tree->Branch("TP_trueY", &trueY);
    tree->Branch("TP_trueZ", &trueZ);
    tree->Branch("TP_signal", &signal);  
  }
  
  void clear() {
    channels.clear(); 
    startT.clear();
    peakT.clear();
    TOT.clear();
    SADC.clear();
    peakADC.clear();
    plane.clear();
    TPC.clear();
    trueX.clear();
    trueY.clear();
    trueZ.clear();
    signal.clear();
  }
};




class TPValTreeWriter : public art::EDAnalyzer {
public:
  explicit TPValTreeWriter(fhicl::ParameterSet const& p);

  TPValTreeWriter(TPValTreeWriter const&) = delete;
  TPValTreeWriter(TPValTreeWriter&&) = delete;
  TPValTreeWriter& operator=(TPValTreeWriter const&) = delete;
  TPValTreeWriter& operator=(TPValTreeWriter&&) = delete;

  
  void analyze(art::Event const& e) override;
  void beginJob() override;
 
private:

  void ResetVariables();
  //function for matching TPs to sim::IDEs 
  std::vector<const sim::IDE*> TPToSimIDEs_Ps(recob::Hit const& hit) const;

  //fcl configurable sampling rate 
  int fADC_SAMPLING_RATE_IN_DTS; // DTS time ticks between adc samples

  // Producer module, configurable from fhicl
  // art::InputTag fRawDigitLabel;
  art::InputTag fTPLabel;
  art::InputTag fGenLabel; //generator label for "signal" particles
  std::string fSimChanLabel; // sim channel label

  bool fSaveMC;
  bool fSaveNeutrino;
  bool fSaveTPs; 
  bool fSaveQs; 

  int foffsetU; //offset for valid TP window for IDE matching (needed for induction sigs).  
  int foffsetV; 
  int foffsetX; 

  //ROOT tree for storing output information 
  TTree* fTree; 
  TTree* fMcTree; 
  TTree* fQTree; 

  // int fEvent, fRun, fSubRun;
  // int fn_TPs; 
  // int fn_noise_TPs_X; 
  // int fn_noise_TPs_U; 
  // int fn_noise_TPs_V; 

  // int fn_signal_TPs_X;
  // int fn_signal_TPs_U;
  // int fn_signal_TPs_V;

  // //total charge for the event based on  #electrons at the readout plane
  // double total_charge_X;
  // double total_charge_U;
  // double total_charge_V; 

  // //charge of IDEs that were detected via hit finding  for the event 
  // double detected_charge_X;
  // double detected_charge_U;
  // double detected_charge_V;


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


  // Monte carlo truth information
  // unsigned int fnParticles; //total number of geant particles 
  // std::vector<int>     fTrackId;
  // std::vector<float>   fMother; 
  // std::vector<float>   fEng;
  // std::vector<float>   fEkin;
  // std::vector<float>   fMass;
  // std::vector<int>     fPdg;
  // std::vector<float>   fP;
  // std::vector<float>   fPx;
  // std::vector<float>   fPy;
  // std::vector<float>   fPz;
  // std::vector<double>  fstartX;
  // std::vector<double>  fstartY;
  // std::vector<double>  fstartZ;
  // std::vector<double>  fendX;
  // std::vector<double>  fendY;
  // std::vector<double>  fendZ;

  // unsigned int fnIDEs; //total number of geant particles 
  // std::vector<uint32_t> fIDE_channel;
  // std::vector<float> fIDE_time;
  // std::vector<int> fIDE_track_id;
  // std::vector<float> fIDE_numElectrons;
  // std::vector<float> fIDE_energy;
  // std::vector<float> fIDE_x;
  // std::vector<float> fIDE_y;
  // std::vector<float> fIDE_z;


  EventDataBuffer fEventData;
  MCTruthDataBuffer fMcTruthData;
  IDEsDataBuffer fIDEsData;
  TPSummaryDataBuffer fTPSummaryData;
  TPsDataBuffer fTPsData;

};

// -- Constructor --
TPValTreeWriter::TPValTreeWriter(fhicl::ParameterSet const& p)
  : EDAnalyzer{p},
  fADC_SAMPLING_RATE_IN_DTS(p.get<int>("ADC_SAMPLING_RATE_IN_DTS",32)),
  fGenLabel(p.get<art::InputTag>("gen_tag")),
  fSimChanLabel(p.get<std::string>("simch_tag", "tpcrawdecoder:simpleSC")),
  // fRawDigitLabel(p.get<std::string>("rawdigit_tag", "tpcrawdecoder:daq")),
  fSaveMC(p.get<bool>("SaveMCInfo",true)), 
  fSaveNeutrino(p.get<bool>("SaveNeutrino",false)),
  fSaveTPs(p.get<bool>("SaveTPInfo",true)),
  fSaveQs(p.get<bool>("SaveQInfo",false)),
  foffsetU(p.get<int>("U_window_end_offset",0)),
  foffsetV(p.get<int>("V_window_end_offset",0)),
  foffsetX(p.get<int>("X_window_end_offset",0))
  {

    if (fSaveTPs) {
      fTPLabel = p.get<art::InputTag>("tp_tag"); // Get TP label from fcl file. 
    }

  }


void TPValTreeWriter::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs; 
  fTree = tfs->make<TTree>("tree", "Trigger Primitives Tree");
  fMcTree = tfs->make<TTree>("mctree", "MC Tree");
  fQTree = tfs->make<TTree>("qtree", "Charge Tree");
  

  fEventData.branch_on(fTree);
  // fTree->Branch("event",&fEvent,"event/i");
  // fTree->Branch("run",&fRun,"run/i");
  // fTree->Branch("subrun",&fSubRun,"subrun/i");

  fTPSummaryData.branch_on(fTree);
  // TPG info
  // fTree->Branch("n_TPs", &fn_TPs, "n_TPs/i");
  // fTree->Branch("n_noise_TPs_X", &fn_noise_TPs_X, "n_noise_TPs_X/i");
  // fTree->Branch("n_noise_TPs_U", &fn_noise_TPs_U, "n_noise_TPs_U/i");
  // fTree->Branch("n_noise_TPs_V", &fn_noise_TPs_V, "n_noise_TPs_V/i");

  // fTree->Branch("n_signal_TPs_X", &fn_signal_TPs_X, "n_signal_TPs_X/i");
  // fTree->Branch("n_signal_TPs_U", &fn_signal_TPs_U, "n_signal_TPs_U/i");
  // fTree->Branch("n_signal_TPs_V", &fn_signal_TPs_V, "n_signal_TPs_V/i");  


  // fTree->Branch("totQ_X", &total_charge_X, "totQ_X/d");
  // fTree->Branch("totQ_U", &total_charge_U, "totQ_U/d");
  // fTree->Branch("totQ_V", &total_charge_V, "totQ_V/d");

  // fTree->Branch("detQ_X", &detected_charge_X, "detQ_X/d");
  // fTree->Branch("detQ_U", &detected_charge_U, "detQ_U/d");
  // fTree->Branch("detQ_V", &detected_charge_V, "detQ_V/d");


  //TP info
  if (fSaveTPs){ 
    fTPsData.branch_on(fTree);
    // fTree->Branch("TP_channel", &fTP_channels);
    // fTree->Branch("TP_startT", &fTP_startT);
    // fTree->Branch("TP_peakT", &fTP_peakT);
    // fTree->Branch("TP_TOT", &fTP_TOT);
    // fTree->Branch("TP_SADC", &fTP_SADC);
    // fTree->Branch("TP_peakADC", &fTP_peakADC);
    // fTree->Branch("TP_plane", &fTP_plane);
    // fTree->Branch("TP_TPC", &fTP_TPC);
    // fTree->Branch("TP_trueX", &fTP_trueX);
    // fTree->Branch("TP_trueY", &fTP_trueY);
    // fTree->Branch("TP_trueZ", &fTP_trueZ);
    // fTree->Branch("TP_signal", &fTP_signal);  
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
  // fTree->Branch("nParticles",&fnParticles,"nParticles/i"); 
  // fTree->Branch("TrackId",&fTrackId);
  // fTree->Branch("Mother",&fMother);
  // fTree->Branch("Pdg",&fPdg);
  // fTree->Branch("Eng",&fEng);
  // fTree->Branch("Ekin",&fEkin);
  // fTree->Branch("Mass",&fMass);
  // fTree->Branch("P",&fP);
  // fTree->Branch("Px",&fPx);
  // fTree->Branch("Py",&fPy);
  // fTree->Branch("Pz",&fPz);
  // fTree->Branch("startX",&fstartX);
  // fTree->Branch("startY",&fstartY);
  // fTree->Branch("startZ",&fstartZ);
  // fTree->Branch("endX",&fendX);
  // fTree->Branch("endY",&fendY);
  // fTree->Branch("endZ",&fendZ);
  fMcTruthData.branch_on(fMcTree);

  //G4 info

  fEventData.branch_on(fMcTree);
  fMcTruthData.branch_on(fMcTree);

  // TODO: add if statement?
  // fMcTree->Branch("event",&fEvent,"event/i");
  // fMcTree->Branch("run",&fRun,"run/i");
  // fMcTree->Branch("subrun",&fSubRun,"subrun/i");

  // fMcTree->Branch("nParticles",&fnParticles,"nParticles/i"); 
  // fMcTree->Branch("TrackId",&fTrackId);
  // fMcTree->Branch("Mother",&fMother);
  // fMcTree->Branch("Pdg",&fPdg);
  // fMcTree->Branch("Eng",&fEng);
  // fMcTree->Branch("Ekin",&fEkin);
  // fMcTree->Branch("Mass",&fMass);
  // fMcTree->Branch("P",&fP);
  // fMcTree->Branch("Px",&fPx);
  // fMcTree->Branch("Py",&fPy);
  // fMcTree->Branch("Pz",&fPz);
  // fMcTree->Branch("startX",&fstartX);
  // fMcTree->Branch("startY",&fstartY);
  // fMcTree->Branch("startZ",&fstartZ);
  // fMcTree->Branch("endX",&fendX);
  // fMcTree->Branch("endY",&fendY);
  // fMcTree->Branch("endZ",&fendZ);  

  fEventData.branch_on(fQTree);
  fIDEsData.branch_on(fQTree);

  // fQTree->Branch("event",&fEvent,"event/i");
  // fQTree->Branch("run",&fRun,"run/i");
  // fQTree->Branch("subrun",&fSubRun,"subrun/i");

  // fQTree->Branch("nIDEs",&fnIDEs,"nIDEs/i");
  // fQTree->Branch("channel",&fIDE_channel);  
  // fQTree->Branch("time",&fIDE_time);  
  // fQTree->Branch("TrackId",&fIDE_track_id);  
  // fQTree->Branch("nElectrons",&fIDE_numElectrons);  
  // fQTree->Branch("energy",&fIDE_energy);  
  // fQTree->Branch("x",&fIDE_x);  
  // fQTree->Branch("y",&fIDE_y);  
  // fQTree->Branch("z",&fIDE_z);  

  std::cout << ">>>>>>>>>>>>>>>>> Saving detector settings" << std::endl;
  auto const& geo = art::ServiceHandle<geo::Geometry>();

  json j;
  j["detector_name"] = geo->DetectorName();

  auto n = tfs->make<TNamed>("settings", j.dump().c_str());
  n->Write();

}

void TPValTreeWriter::analyze(art::Event const& e) {
  ResetVariables();  // initialise/reset all variables

  // fRun = e.run();
  // fSubRun = e.subRun();
  // fEvent = e.id().event();

  fEventData.event = e.id().event();
  fEventData.run = e.run();
  fEventData.subrun = e.subRun();


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
          // total_charge_X += ide.numElectrons;
          fTPSummaryData.total_charge_X += ide.numElectrons;
        }
        else if (plane == geo::kU) {
          // total_charge_U += ide.numElectrons;
          fTPSummaryData.total_charge_U += ide.numElectrons;
        }
        else if (plane == geo::kV) {
          // total_charge_V += ide.numElectrons;
          fTPSummaryData.total_charge_V += ide.numElectrons;
        }
      }
    }
  }

  if (fSaveTPs){

    // Process the Trigger Primitives
    auto tpHandle = e.getValidHandle<std::vector<dunedaq::trgdataformats::TriggerPrimitive>>(fTPLabel);
    // auto prov = tpHandle.provenance();
    auto const& tps = *tpHandle;
    // fn_TPs = tps.size();
    fTPSummaryData.num_TPs = tps.size();

    for (const auto& tp : tps) {
      auto plane =  wireReadout.ROPtoWirePlanes(wireReadout.ChannelToROP(tp.channel)).at(0).Plane;

      // Create recob::Hit for TP (with a larger time window for induction planes)
      float tp_rms = (plane == 2) ? (tp.time_over_threshold / this->fADC_SAMPLING_RATE_IN_DTS) * 2 : (tp.time_over_threshold / this->fADC_SAMPLING_RATE_IN_DTS) * 4;
      recob::Hit ThisHit(
                        static_cast<raw::ChannelID_t>(tp.channel),        //channel
                        static_cast<raw::TDCtick_t>(tp.time_start / this->fADC_SAMPLING_RATE_IN_DTS), //start tick
                        static_cast<raw::TDCtick_t>((tp.time_start / this->fADC_SAMPLING_RATE_IN_DTS) + (tp.time_over_threshold / this->fADC_SAMPLING_RATE_IN_DTS)), //end tick
                        static_cast<float>((tp.time_peak / this->fADC_SAMPLING_RATE_IN_DTS)), // peak time 
                        static_cast<float>((tp.time_over_threshold / this->fADC_SAMPLING_RATE_IN_DTS) * 0.5), // sigma peak time 
                        static_cast<float>(tp_rms), // rms 
                        static_cast<float>(tp.adc_peak), //peak amplitude 
                        static_cast<float>(0), // sigma peak amplitude 
                        static_cast<float>(tp.adc_integral), //ROI SADC 
                        static_cast<float>(tp.adc_integral), //hit SADC 
                        static_cast<float>(0), //hit integral 
                        static_cast<float>(0), //hit sigma integral 
                        static_cast<short int>(0), //multiplicity 
                        static_cast<short int>(0), //local index 
                        static_cast<float>(0), //goodness of fit 
                        static_cast<int>(0), //degrees of freedom 
                        geo::View_t(plane), //view 
                        wireReadout.SignalType(plane), //sig type 
                        wireReadout.ChannelToWire(tp.channel).front() //wire ID 
                        );

      // Avoid reprocessing the same IDEs for different windows
      std::unordered_set<const sim::IDE*> processedIDEs;

      std::vector<const sim::IDE*> ides;
      try { 
        //use standard backtracker based on peak time for RS/ST & custom one based on end/start times for AbsRS since the peak time is offset 
        //it actually gives worse performance for RS & ST so just use the same approach for everything 
        //ides = fAbsRS ? TPToSimIDEs_Ps(ThisHit) : bt_serv->HitToSimIDEs_Ps(clockData,ThisHit);
        ides = TPToSimIDEs_Ps(ThisHit); 
      } catch(...) {}

      
      // If there are associated IDEs, save the true coordinates
      if (!ides.empty()) {

        switch (plane){
          case geo::kW:
            fTPSummaryData.num_signal_TPs_X++;
            break;
          case geo::kU:
            fTPSummaryData.num_signal_TPs_U++;
            break;
          case geo::kV:
            fTPSummaryData.num_signal_TPs_V++;
            break;
          default:
            break;
        }

        // if (plane == geo::kW) fn_signal_TPs_X++;
        // else if (plane == geo::kU) fn_signal_TPs_U++;
        // else if (plane == geo::kV) fn_signal_TPs_V++;

        // Loop over detected IDEs and avoid double-counting using a set
        for (const sim::IDE* ide_ptr : ides) {
          if (processedIDEs.find(ide_ptr) == processedIDEs.end()) {
            processedIDEs.insert(ide_ptr);
            float ide_numElectrons = ide_ptr->numElectrons;

            if (plane == geo::kW) {
              // detected_charge_X += ide_numElectrons;
              fTPSummaryData.detected_charge_X += ide_numElectrons;
            } else if (plane == geo::kU) {
              // detected_charge_U += ide_numElectrons;
              fTPSummaryData.detected_charge_U += ide_numElectrons;
            } else if (plane == geo::kV) {
              // detected_charge_V += ide_numElectrons;
              fTPSummaryData.detected_charge_V += ide_numElectrons;
            }
          }
        }
        std::vector<double> xyz = bt_serv->SimIDEsToXYZ(ides);
        fTPsData.trueX.push_back(xyz[0]);
        fTPsData.trueY.push_back(xyz[1]);
        fTPsData.trueZ.push_back(xyz[2]);
        fTPsData.signal.push_back(1);
      } else {

        // If no IDEs, count as noise TP
        switch (plane){
          case geo::kW:
            fTPSummaryData.num_noise_TPs_X++;
            break;
          case geo::kU:
            fTPSummaryData.num_noise_TPs_U++;
            break;
          case geo::kV:
            fTPSummaryData.num_noise_TPs_V++;
            break;
          default:
            break;
        }
        // if (plane == geo::kW) fTPSummaryData.n_noise_TPs_X++;
        // else if (plane == geo::kU) fTPSummaryData.n_noise_TPs_U++;
        // else if (plane == geo::kV) fTPSummaryData.n_noise_TPs_V++;
        
        fTPsData.signal.push_back(0);
        fTPsData.trueX.push_back(-1);
        fTPsData.trueY.push_back(-1);
        fTPsData.trueZ.push_back(-1);
        
      }
      // Store TP info 
      fTPsData.channels.push_back(tp.channel);
      fTPsData.startT.push_back(tp.time_start / this->fADC_SAMPLING_RATE_IN_DTS);
      fTPsData.peakT.push_back((tp.time_peak / this->fADC_SAMPLING_RATE_IN_DTS));
      fTPsData.TOT.push_back(tp.time_over_threshold / this->fADC_SAMPLING_RATE_IN_DTS);
      fTPsData.SADC.push_back(tp.adc_integral);
      fTPsData.peakADC.push_back(tp.adc_peak);
      fTPsData.plane.push_back(plane);
      fTPsData.TPC.push_back(wireReadout.ROPtoTPCs(wireReadout.ChannelToROP(tp.channel)).at(0).TPC);
    }
  }


  //access MC particles through truth label instead of G4 handle, 
  //only interested in "signal" particles &  don't want to iterate over & store info for the avalanche of radiologicals
    
  if (fSaveMC){
    auto mc_truth =  e.getValidHandle< std::vector<simb::MCTruth> >(fGenLabel);
 
    if (mc_truth.isValid()){
      for (const auto& truth : *mc_truth) {

        // fnParticles = truth.NParticles();    
        fMcTruthData.num_particles = truth.NParticles();    

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
          fMcTruthData.TrackId.push_back( trueParticle.TrackId());
          fMcTruthData.Mother.push_back( trueParticle.Mother());
          fMcTruthData.Eng.push_back( trueParticle.E() );
          fMcTruthData.Pdg.push_back( trueParticle.PdgCode());
          fMcTruthData.Ekin.push_back( Ekin ); 
          fMcTruthData.Mass.push_back( trueParticle.Mass() );
          fMcTruthData.P.push_back( trueParticle.P()) ;
          fMcTruthData.Px.push_back( trueParticle.Px());
          fMcTruthData.Py.push_back( trueParticle.Py());
          fMcTruthData.Pz.push_back( trueParticle.Pz());
          fMcTruthData.startX.push_back( trueParticle.Vx());
          fMcTruthData.startY.push_back( trueParticle.Vy());
          fMcTruthData.startZ.push_back( trueParticle.Vz());
          fMcTruthData.endX.push_back( trueParticle.EndPosition()[0]);
          fMcTruthData.endY.push_back( trueParticle.EndPosition()[1]);
          fMcTruthData.endZ.push_back( trueParticle.EndPosition()[2]);
        }
      }
    }
  }

  // Process the SimChannels just once
  for (const sim::SimChannel &sc : *simchannels) {
    auto tdcidemap = sc.TDCIDEMap();

    // Sum charge for each plane
    for (const auto& [t, tdcide] : tdcidemap) {
      for (const auto& ide : tdcide) {
          fIDEsData.channel.push_back(sc.Channel());
          fIDEsData.time.push_back(t);
          fIDEsData.track_id.push_back(ide.trackID);
          fIDEsData.numElectrons.push_back(ide.numElectrons);
          fIDEsData.energy.push_back(ide.energy);
          fIDEsData.x.push_back(ide.x);
          fIDEsData.y.push_back(ide.y);
          fIDEsData.z.push_back(ide.z);
      }
      fIDEsData.entries += tdcide.size();
    }


  }


  fTree->Fill();
  fMcTree->Fill();
  fQTree->Fill();

}


//TP - ide matching based on hit start and end times rather than peak ADC time 
std::vector<const sim::IDE*> TPValTreeWriter::TPToSimIDEs_Ps(recob::Hit const& hit) const
{
  std::vector<const sim::IDE*> retVec;
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  
  //IDEs are computed based on the time during which the ionization cloud hits the readout plane. For ind. signals this corresponds to the inflection point
  //need to introduce an offset for this effect to be properly accounted for 
  int offset = 0;

  // if (hit.View() == geo::kU) offset = foffsetU; //fcl configurable offset for induction TPs 
  // else if (hit.View() == geo::kV) offset = foffsetV;  

  switch (hit.View()) {
    case geo::kU:
      offset = foffsetU;
      break;
    case geo::kV:
      offset = foffsetV;
      break;
    case geo::kW:
      offset = foffsetX;
      break;
    default:
      break;
  }


  // Adjust the start and end time based on the hit, this may need to be expanded or modified
  int start_tdc = hit.StartTick() + offset; //clockData.TPCTick2TDC(hit.PeakTimeMinusRMS(fHitTimeRMS));
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



void TPValTreeWriter::ResetVariables()
{
  // fEvent = fRun = fSubRun = -1; 
  // fn_TPs = -1; 


  // total_charge_X = total_charge_U = total_charge_V = 0; 


  // detected_charge_X = detected_charge_U = detected_charge_V = 0; 

  // //initialise number of noise and signal TPs for this event 
  // fn_noise_TPs_X = fn_noise_TPs_U = fn_noise_TPs_V =  0;
  // fn_signal_TPs_X = fn_signal_TPs_U = fn_signal_TPs_V =  0;

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


  // // GEANT4 info
  // fnParticles = 0;
  // fTrackId.clear();
  // fMother.clear();
  // fEng.clear();
  // fEkin.clear();
  // fMass.clear();
  // fPdg.clear();
  // fP.clear();
  // fPx.clear();
  // fPy.clear();
  // fPz.clear();
  // fstartX.clear();
  // fstartY.clear();
  // fstartZ.clear();
  // fendX.clear();
  // fendY.clear();
  // fendZ.clear();


  // fnIDEs = 0;
  // fIDE_channel.clear();
  // fIDE_time.clear();
  // fIDE_track_id.clear();
  // fIDE_numElectrons.clear();
  // fIDE_energy.clear();
  // fIDE_x.clear();
  // fIDE_y.clear();
  // fIDE_z.clear();


  fEventData.clear();
  fMcTruthData.clear();
  fIDEsData.clear();
  fTPSummaryData.clear();
  fTPsData.clear();
}

}

DEFINE_ART_MODULE(duneana::TPValTreeWriter)
