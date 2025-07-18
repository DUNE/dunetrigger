////////////////////////////////////////////////////////////////////////
// Class:       TriggerAnaTree
// Plugin Type: analyzer (Unknown Unknown)
// File:        TriggerAnaTree_module.cc
//
// Generated at Fri Aug 30 14:50:19 2024 by jierans using cetskelgen
// from cetlib version 3.18.02.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "detdataformats/trigger/TriggerActivityData.hpp"
#include "detdataformats/trigger/TriggerCandidateData.hpp"
#include "detdataformats/trigger/TriggerPrimitive.hpp"

#include "lardataobj/Simulation/SimChannel.h"
#include "larcore/Geometry/WireReadout.h"
#include "larcoreobj/SimpleTypesAndConstants/readout_types.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include <TDirectory.h>
#include <TFile.h>
#include <TTree.h>

using dunedaq::trgdataformats::TriggerActivityData;
using dunedaq::trgdataformats::TriggerCandidateData;
using dunedaq::trgdataformats::TriggerPrimitive;

namespace dunetrigger{
  class TriggerAnaTree;
}

class dunetrigger::TriggerAnaTree : public art::EDAnalyzer {
public:
  typedef struct ChannelInfo {
    unsigned int rop_id;
    int view;
    unsigned int tpcset_id;
  } matching_info;

  explicit TriggerAnaTree(fhicl::ParameterSet const &p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  TriggerAnaTree(TriggerAnaTree const &) = delete;
  TriggerAnaTree(TriggerAnaTree &&) = delete;
  TriggerAnaTree &operator=(TriggerAnaTree const &) = delete;
  TriggerAnaTree &operator=(TriggerAnaTree &&) = delete;

  // Required functions.
  void beginJob() override;
  void analyze(art::Event const &e) override;
  // void endJob() override;

private:
  art::ServiceHandle<art::TFileService> tfs;
  std::map<std::string, TTree *> tree_map;
  // buffers for writing to ROOT Trees
  unsigned int fEventID;
  int fRun;
  int fSubRun;
  int fAssnIdx;

  std::unordered_map<int, int> trkId_to_truthBlockId;
  std::map<std::string, dunedaq::trgdataformats::TriggerPrimitive> tp_bufs;
  std::map<std::string, ChannelInfo> tp_channel_info_bufs;
  std::map<std::string, dunedaq::trgdataformats::TriggerActivityData> ta_bufs;
  std::map<std::string, dunedaq::trgdataformats::TriggerCandidateData> tc_bufs;

  bool dump_tp, dump_ta, dump_tc;

  void make_tp_tree_if_needed(std::string tag, bool assn = false);
  void make_ta_tree_if_needed(std::string tag, bool assn = false);
  void make_tc_tree_if_needed(std::string tag);

  ChannelInfo get_channel_info_for_channel(geo::WireReadoutGeom const *geom, int channel);

  bool dump_mctruths;
  TTree *mctruth_tree;
  int mctruth_pdg;
  std::string mctruth_process;
  int mctruth_status, mctruth_id;
  std::string mctruth_gen_name;
  double mctruth_x, mctruth_y, mctruth_z;
  double mctruth_Px, mctruth_Py, mctruth_Pz;
  double mctruth_en;

  bool dump_mcparticles;
  TTree *mcparticle_tree;
  int mcparticle_pdg;
  std::string mcparticle_process;
  int mcparticle_status, mcparticle_trackid, mcparticle_truthid;
  std::string mcparticle_gen_name;
  double mcparticle_x, mcparticle_y, mcparticle_z;
  double mcparticle_Px, mcparticle_Py, mcparticle_Pz;
  double mcparticle_en;

  bool dump_simides;
  std::string simchannel_tag;
  TTree *simide_tree;
  unsigned int sim_channel_id;
  unsigned short tdc;
  float ide_numElectrons, ide_energy, ide_x, ide_y, ide_z;
  int ide_trkId, ide_origTrkId;
};

dunetrigger::TriggerAnaTree::TriggerAnaTree(fhicl::ParameterSet const &p)
    : EDAnalyzer{p}, dump_tp(p.get<bool>("dump_tp")),
      dump_ta(p.get<bool>("dump_ta")), dump_tc(p.get<bool>("dump_tc")),
      dump_mctruths(p.get<bool>("dump_mctruths", true)),
      dump_mcparticles(p.get<bool>("dump_mcparticles", true)),
      dump_simides(p.get<bool>("dump_simides", true)),
      simchannel_tag(
          p.get<std::string>("simchannel_tag", "tpcrawdecoder:simpleSC"))
// More initializers here.
{}

void dunetrigger::TriggerAnaTree::beginJob() {
  if (dump_mctruths) {
    mctruth_tree = tfs->make<TTree>("mctruths", "mctruths");
    mctruth_tree->Branch("Event", &fEventID, "Event/i");
    mctruth_tree->Branch("Run", &fRun, "Run/i");
    mctruth_tree->Branch("SubRun", &fSubRun, "SubRun/i");
    mctruth_tree->Branch("block_id", &mctruth_id);
    mctruth_tree->Branch("pdg", &mctruth_pdg);
    mctruth_tree->Branch("generator_name", &mctruth_gen_name);
    mctruth_tree->Branch("status_code", &mctruth_status);
    mctruth_tree->Branch("x", &mctruth_x);
    mctruth_tree->Branch("y", &mctruth_y);
    mctruth_tree->Branch("z", &mctruth_z);
    mctruth_tree->Branch("Px", &mctruth_Px);
    mctruth_tree->Branch("Py", &mctruth_Py);
    mctruth_tree->Branch("Pz", &mctruth_Pz);
    mctruth_tree->Branch("en", &mctruth_en);
    mctruth_tree->Branch("process", &mctruth_process);
  }
  if (dump_mcparticles) {
    mcparticle_tree = tfs->make<TTree>("mcparticles", "mcparticles");
    mcparticle_tree->Branch("Event", &fEventID, "Event/i");
    mcparticle_tree->Branch("Run", &fRun, "Run/i");
    mcparticle_tree->Branch("SubRun", &fSubRun, "SubRun/i");
    mcparticle_tree->Branch("pdg", &mcparticle_pdg);
    mcparticle_tree->Branch("generator_name", &mcparticle_gen_name);
    mcparticle_tree->Branch("status_code", &mcparticle_status);
    mcparticle_tree->Branch("track_id", &mcparticle_trackid);
    mcparticle_tree->Branch("truth_id", &mcparticle_truthid);
    mcparticle_tree->Branch("x", &mcparticle_x);
    mcparticle_tree->Branch("y", &mcparticle_y);
    mcparticle_tree->Branch("z", &mcparticle_z);
    mcparticle_tree->Branch("Px", &mcparticle_Px);
    mcparticle_tree->Branch("Py", &mcparticle_Py);
    mcparticle_tree->Branch("Pz", &mcparticle_Pz);
    mcparticle_tree->Branch("en", &mcparticle_en);
    mcparticle_tree->Branch("process", &mcparticle_process);
  }
  if (dump_simides) {
    simide_tree = tfs->make<TTree>("simides", "simides");
    simide_tree->Branch("Event", &fEventID, "Event/i");
    simide_tree->Branch("Run", &fRun, "Run/i");
    simide_tree->Branch("SubRun", &fSubRun, "SubRun/i");
    simide_tree->Branch("ChannelID", &sim_channel_id, "ChannelID/i");
    simide_tree->Branch("Timestamp", &tdc, "Timestamp/s");
    simide_tree->Branch("numElectrons", &ide_numElectrons, "ide_numElectrons/i");
    simide_tree->Branch("Energy", &ide_energy, "ide_energy/F");
    simide_tree->Branch("x", &ide_x, "ide_x/F");
    simide_tree->Branch("y", &ide_y, "ide_y/F");
    simide_tree->Branch("z", &ide_z, "ide_z/F");
    simide_tree->Branch("trackID", &ide_trkId, "ide_trkId/I");
    simide_tree->Branch("origTrackID", &ide_origTrkId, "ide_origTrkId/I");
  }
}

void dunetrigger::TriggerAnaTree::analyze(art::Event const &e) {
  fRun = e.run();
  fSubRun = e.subRun();
  fEventID = e.id().event();
  // get a service handle for geometry
  geo::WireReadoutGeom const *geom =
      &art::ServiceHandle<geo::WireReadout>()->Get();
  if (dump_mctruths) {
    std::vector<art::Handle<std::vector<simb::MCTruth>>> mctruthHandles =
        e.getMany<std::vector<simb::MCTruth>>();
    int truth_block_counter = 0;
    trkId_to_truthBlockId.clear();
    for (auto const &mctruthHandle : mctruthHandles) {
      std::string generator_name =
          mctruthHandle.provenance()->inputTag().label();
      // NOTE: here we are making an assumption that the geant4 stage's process name is largeant.
      // This should be safe mostly.
      art::FindManyP<simb::MCParticle> assns(mctruthHandle, e, "largeant");
      for (size_t i = 0; i < mctruthHandle->size(); i++){
        const simb::MCTruth &truthblock = *art::Ptr<simb::MCTruth>(mctruthHandle, i);
        int nparticles = truthblock.NParticles();
        if (dump_mcparticles){
          std::vector<art::Ptr<simb::MCParticle>> matched_mcparts = assns.at(i);
          for (art::Ptr<simb::MCParticle> mcpart: matched_mcparts) {
            trkId_to_truthBlockId[mcpart->TrackId()] = truth_block_counter;
          }
        }
        for (int ipart = 0; ipart < nparticles; ipart++) {
          const simb::MCParticle &part = truthblock.GetParticle(ipart);
          mctruth_id = truth_block_counter;
          mctruth_pdg = part.PdgCode();
          mctruth_gen_name = generator_name;
          mctruth_status = part.StatusCode();
          mctruth_process = part.Process();
          mctruth_x = part.Vx();
          mctruth_y = part.Vy();
          mctruth_z = part.Vz();
          mctruth_Px = part.Px();
          mctruth_Py = part.Py();
          mctruth_Pz = part.Pz();
          mctruth_en = part.E();
          mctruth_tree->Fill();
        }
        truth_block_counter++;
      }
    }
  }
  if (dump_mcparticles) {
    std::vector<art::Handle<std::vector<simb::MCParticle>>> mcparticleHandles =
        e.getMany<std::vector<simb::MCParticle>>();
    for (auto const &mcparticleHandle : mcparticleHandles) {
      std::string generator_name =
          mcparticleHandle.provenance()->inputTag().label();
      for (const simb::MCParticle &part : *mcparticleHandle) {
        mcparticle_pdg = part.PdgCode();
        mcparticle_gen_name = generator_name;
        mcparticle_status = part.StatusCode();
        mcparticle_trackid = part.TrackId();
        mcparticle_truthid = dump_mctruths ? trkId_to_truthBlockId.at(part.TrackId()) : -1;
        mcparticle_process = part.Process();
        mcparticle_x = part.Vx();
        mcparticle_y = part.Vy();
        mcparticle_z = part.Vz();
        mcparticle_Px = part.Px();
        mcparticle_Py = part.Py();
        mcparticle_Pz = part.Pz();
        mcparticle_en = part.E();
        mcparticle_tree->Fill();
      }
    }
  }
  if (dump_simides) {

    auto simchannels =
        e.getValidHandle<std::vector<sim::SimChannel>>(simchannel_tag);
    for (const sim::SimChannel &sc : *simchannels) {
      sim_channel_id = sc.Channel();
      sim::SimChannel::TDCIDEs_t const &tdcidemap = sc.TDCIDEMap();
      for (const sim::TDCIDE &tdcide : tdcidemap) {
        tdc = tdcide.first;
        for (sim::IDE ide : tdcide.second) {
          ide_numElectrons = ide.numElectrons;
          ide_energy = ide.energy;
          ide_x = ide.x;
          ide_y = ide.y;
          ide_z = ide.z;
          ide_trkId = ide.trackID;
          ide_origTrkId = ide.origTrackID;
          simide_tree->Fill();
        }
      }
    }
  }
  if (dump_tp) {
    std::vector<art::Handle<std::vector<TriggerPrimitive>>> tpHandles =
        e.getMany<std::vector<TriggerPrimitive>>();
    for (auto const &tpHandle : tpHandles) {
      std::string tag = tpHandle.provenance()->inputTag().encode();
      std::string map_tag = "tp/" + tag;
      make_tp_tree_if_needed(tag);
      for (const TriggerPrimitive &tp : *tpHandle) {
        tp_bufs[map_tag] = tp;
        tp_channel_info_bufs[map_tag] = get_channel_info_for_channel(geom, tp.channel);
        tree_map[map_tag]->Fill();
      }
    }
  }

  if (dump_ta) {
    std::vector<art::Handle<std::vector<TriggerActivityData>>> taHandles =
        e.getMany<std::vector<TriggerActivityData>>();
    for (auto const &taHandle : taHandles) {
      art::FindManyP<TriggerPrimitive> assns(
          taHandle, e, taHandle.provenance()->moduleLabel());
      std::string tag = taHandle.provenance()->inputTag().encode();
      std::string map_tag = "ta/" + tag;
      make_ta_tree_if_needed(tag);
      for (unsigned int i = 0; i < taHandle->size(); i++) {
        const TriggerActivityData &ta =
            *art::Ptr<TriggerActivityData>(taHandle, i);
        if (assns.isValid()) {
          art::InputTag ta_input_tag = taHandle.provenance()->inputTag();
          std::string tpInTaTag =
              art::InputTag(ta_input_tag.label(),
                            ta_input_tag.instance() + "inTAs",
                            ta_input_tag.process())
                  .encode();
          std::string map_tpInTaTag = "tp/" + tpInTaTag;
          make_tp_tree_if_needed(tpInTaTag, true);
          fAssnIdx = i;
          std::vector<art::Ptr<TriggerPrimitive>> matched_tps = assns.at(i);
          for (art::Ptr<TriggerPrimitive> tp : matched_tps) {
            tp_bufs[map_tpInTaTag] = *tp;
            tp_channel_info_bufs[map_tag] = get_channel_info_for_channel(geom, tp->channel);
            tree_map[map_tpInTaTag]->Fill();
          }
        }
        ta_bufs[map_tag] = ta;
        tree_map[map_tag]->Fill();
      }
    }
  }

  if (dump_tc) {
    std::vector<art::Handle<std::vector<TriggerCandidateData>>> tcHandles =
        e.getMany<std::vector<TriggerCandidateData>>();
    for (auto const &tcHandle : tcHandles) {
      art::FindManyP<TriggerActivityData> assns(
          tcHandle, e, tcHandle.provenance()->moduleLabel());
      std::string tag = tcHandle.provenance()->inputTag().encode();
      std::string map_tag = "tc/" + tag;
      make_tc_tree_if_needed(tag);
      for (unsigned int i = 0; i < tcHandle->size(); i++) {
        const TriggerCandidateData &tc =
            *art::Ptr<TriggerCandidateData>(tcHandle, i);
        if (assns.isValid()) {
          art::InputTag tc_input_tag = tcHandle.provenance()->inputTag();
          std::string taInTcTag =
              art::InputTag(tc_input_tag.label(),
                            tc_input_tag.instance() + "inTCs",
                            tc_input_tag.process())
                  .encode();
          std::string map_taInTcTag = "ta/" + taInTcTag;
          make_ta_tree_if_needed(taInTcTag, true);
          fAssnIdx = i;
          std::vector<art::Ptr<TriggerActivityData>> matched_tas = assns.at(i);
          for (art::Ptr<TriggerActivityData> ta : matched_tas) {
            ta_bufs[map_taInTcTag] = *ta;
            tree_map[map_taInTcTag]->Fill();
          }
        }
        tc_bufs[map_tag] = tc;
        tree_map[map_tag]->Fill();
      }
    }
  }
}

void dunetrigger::TriggerAnaTree::make_tp_tree_if_needed(std::string tag, bool assn) {
  std::string map_tag = "tp/" + tag;
  if (!tree_map.count(map_tag)) {
    art::TFileDirectory tp_dir =
        tfs->mkdir("TriggerPrimitives", "Trigger Primitive Trees");
    std::cout << "Creating new TTree for " << tag << std::endl;

    // Replace ":" with "_" in TTree names so that they can be used in ROOT's
    // intepreter
    std::string tree_name = tag;
    std::replace(tree_name.begin(), tree_name.end(), ':', '_');
    TTree *tree = tp_dir.make<TTree>(tree_name.c_str(), tree_name.c_str());
    tree_map[map_tag] = tree;
    TriggerPrimitive &tp = tp_bufs[map_tag];
    tree->Branch("Event", &fEventID, "Event/i");
    tree->Branch("Run", &fRun, "Run/i");
    tree->Branch("SubRun", &fSubRun, "SubRun/i");
    tree->Branch("version", &tp.version);
    tree->Branch("time_start", &tp.time_start);
    tree->Branch("time_peak", &tp.time_peak);
    tree->Branch("time_over_threshold", &tp.time_over_threshold);
    tree->Branch("channel", &tp.channel);
    tree->Branch("adc_integral", &tp.adc_integral);
    tree->Branch("adc_peak", &tp.adc_peak);
    tree->Branch("detid", &tp.detid);
    tree->Branch("type", &tp.type, "type/I");
    tree->Branch("algorithm", &tp.algorithm, "algorithm/I");
    ChannelInfo &chinfo = tp_channel_info_bufs[map_tag];
    tree->Branch("ropid", &chinfo.rop_id, "ropid/i");
    tree->Branch("view", &chinfo.view, "view/I");
    tree->Branch("TPCSetID", &chinfo.tpcset_id, "TPCSetID/i");
    if (assn)
      tree->Branch("TAnumber", &fAssnIdx, "TAnumber/I");
  }
}

void dunetrigger::TriggerAnaTree::make_ta_tree_if_needed(std::string tag, bool assn) {
  std::string map_tag = "ta/" + tag;
  if (!tree_map.count(map_tag)) {
    art::TFileDirectory ta_dir =
        tfs->mkdir("TriggerActivities", "Trigger Activity Trees");
    std::cout << "Creating new TTree for " << tag << std::endl;
    // Replace ":" with "_" in TTree names so that they can be used in ROOT's
    // intepreter
    std::string tree_name = tag;
    std::replace(tree_name.begin(), tree_name.end(), ':', '_');
    TTree *tree = ta_dir.make<TTree>(tree_name.c_str(), tree_name.c_str());
    tree_map[map_tag] = tree;
    TriggerActivityData &ta = ta_bufs[map_tag];
    tree->Branch("Event", &fEventID, "Event/i");
    tree->Branch("Run", &fRun, "Run/i");
    tree->Branch("SubRun", &fSubRun, "SubRun/i");
    tree->Branch("version", &ta.version);
    tree->Branch("time_start", &ta.time_start);
    tree->Branch("time_end", &ta.time_end);
    tree->Branch("time_peak", &ta.time_peak);
    tree->Branch("time_activity", &ta.time_activity);
    tree->Branch("channel_start", &ta.channel_start);
    tree->Branch("channel_end", &ta.channel_end);
    tree->Branch("channel_peak", &ta.channel_peak);
    tree->Branch("adc_integral", &ta.adc_integral);
    tree->Branch("adc_peak", &ta.adc_peak);
    tree->Branch("detid", &ta.detid);
    // HACK: assuming enums are ints here
    tree->Branch("type", &ta.type, "type/I");
    tree->Branch("algorithm", &ta.algorithm, "algorithm/I");
    if (assn)
      tree->Branch("TCnumber", &fAssnIdx, "TCnumber/I");
  }
}

void dunetrigger::TriggerAnaTree::make_tc_tree_if_needed(std::string tag) {
  std::string map_tag = "tc/" + tag;
  if (!tree_map.count(map_tag)) {
    art::TFileDirectory tc_dir =
        tfs->mkdir("TriggerCandidates", "Trigger Candidate Trees");
    std::cout << "Creating new TTree for " << tag << std::endl;
    // Replace ":" with "_" in TTree names so that they can be used in ROOT's
    // intepreter
    std::string tree_name = tag;
    std::replace(tree_name.begin(), tree_name.end(), ':', '_');
    TTree *tree = tc_dir.make<TTree>(tree_name.c_str(), tree_name.c_str());
    tree_map[map_tag] = tree;
    TriggerCandidateData &tc = tc_bufs[map_tag];
    tree->Branch("Event", &fEventID, "Event/i");
    tree->Branch("Run", &fRun, "Run/i");
    tree->Branch("SubRun", &fSubRun, "SubRun/i");
    tree->Branch("version", &tc.version);
    tree->Branch("time_start", &tc.time_start);
    tree->Branch("time_end", &tc.time_end);
    tree->Branch("time_candidate", &tc.time_candidate);
    tree->Branch("detid", &tc.detid);
    tree->Branch("type", &tc.type, "type/I");
    tree->Branch("algorithm", &tc.algorithm, "algorithm/I");
  }
}

dunetrigger::TriggerAnaTree::ChannelInfo dunetrigger::TriggerAnaTree::get_channel_info_for_channel(geo::WireReadoutGeom const *geom, int channel) {
  readout::ROPID rop = geom->ChannelToROP(channel);
  ChannelInfo result;
  result.rop_id = rop.ROP;
  result.tpcset_id = rop.asTPCsetID().TPCset;
  result.view = geom->View(rop);
  return result;
}

DEFINE_ART_MODULE(dunetrigger::TriggerAnaTree)
