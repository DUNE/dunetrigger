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
#include "dunetrigger/TriggerSim/TPAlgTools/TPAlgTPCTool.hh"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "detdataformats/trigger/TriggerActivityData.hpp"
#include "detdataformats/trigger/TriggerCandidateData.hpp"
#include "detdataformats/trigger/TriggerPrimitive.hpp"

#include "larcore/Geometry/WireReadout.h"
#include "larcoreobj/SimpleTypesAndConstants/readout_types.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larsim/MCCheater/BackTrackerService.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include <TDirectory.h>
#include <TFile.h>
#include <TTree.h>

using dunedaq::trgdataformats::TriggerActivityData;
using dunedaq::trgdataformats::TriggerCandidateData;
using dunedaq::trgdataformats::TriggerPrimitive;

#define INVALID -99999

namespace dunetrigger {
class TriggerAnaTree;

struct ChannelInfo {
  unsigned int rop_id;
  int view;
  unsigned int tpcset_id;
};

// Adapts similar storage optimization techniques as TPv2, but expanding all
// bitfields.
struct TriggerPrimitiveBuffer {
  uint8_t version;
  uint8_t flag;
  uint8_t detid;

  uint32_t channel;
  uint16_t samples_over_threshold;
  uint64_t time_start;
  uint16_t samples_to_peak;
  uint32_t adc_integral;
  uint16_t adc_peak;

  int bt_primary_track_id;
  double bt_primary_track_numelectron_frac, bt_primary_track_energy_frac;
  double bt_edep;
  double bt_numelectrons;
  double bt_x, bt_y, bt_z;
  double bt_primary_x, bt_primary_y, bt_primary_z;

  void from_tp(const TriggerPrimitive &tp);
  void populate_backtracking_info(const std::vector<sim::IDE> &ides);
};

} // namespace dunetrigger

class dunetrigger::TriggerAnaTree : public art::EDAnalyzer {
public:
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
  uint32_t fEventID;
  uint32_t fRun;
  uint32_t fSubRun;
  size_t fAssnIdx;

  std::unordered_map<int, int> trkId_to_truthBlockId;
  std::map<std::string, TriggerPrimitiveBuffer> tp_bufs;
  std::map<std::string, ChannelInfo> tp_channel_info_bufs;
  std::map<std::string, TriggerActivityData> ta_bufs;
  std::map<std::string, TriggerCandidateData> tc_bufs;
  std::map<int, double> track_en_sums;
  std::map<int, double> track_electron_sums;

  bool dump_tp, dump_ta, dump_tc;

  bool tp_backtracking;
  int fOffsetU, fOffsetV, fOffsetX;

  void make_tp_tree_if_needed(std::string tag, bool assn = false);
  void make_ta_tree_if_needed(std::string tag, bool assn = false);
  void make_tc_tree_if_needed(std::string tag);

  std::vector<sim::IDE> match_simides_to_tps(const TriggerPrimitiveBuffer &tp,
                                             const ChannelInfo &chinfo) const;

  ChannelInfo get_channel_info_for_channel(geo::WireReadoutGeom const *geom,
                                           int channel);

  bool dump_mctruths;
  TTree *mctruth_tree;
  int mctruth_pdg;
  std::string mctruth_process;
  int mctruth_status, mctruth_id, mctruth_trackid;
  std::string mctruth_gen_name;
  double mctruth_x, mctruth_y, mctruth_z;
  double mctruth_Px, mctruth_Py, mctruth_Pz;
  double mctruth_en;
  TTree *mcneutrino_tree;
  int mcneutrino_nupdg, mcneutrino_leptonpdg, mcneutrino_ccnc, mcneutrino_mode,
      mcneutrino_iteractionType, mcneutrino_target, mcneutrino_hitnuc,
      mcneutrino_hitquark;
  double mcneutrino_w, mcneutrino_x, mcneutrino_y, mcneutrino_qsqr,
      mcneutrino_pt, mcneutrino_theta;

  bool dump_mcparticles;
  TTree *mcparticle_tree;
  int mcparticle_pdg;
  std::string mcparticle_process;
  int mcparticle_status, mcparticle_trackid, mcparticle_truthid,
      mcparticle_mother;
  std::string mcparticle_gen_name;
  double mcparticle_x, mcparticle_y, mcparticle_z, mcparticle_t;
  double mcparticle_end_x, mcparticle_end_y, mcparticle_end_z, mcparticle_end_t;
  double mcparticle_Px, mcparticle_Py, mcparticle_Pz;
  double mcparticle_en;
  double mcparticle_edep, mcparticle_numelectrons;
  double mcparticle_shower_edep, mcparticle_shower_numelectrons;

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
      tp_backtracking(p.get<bool>("tp_backtracking", false)),
      fOffsetU(p.get<int>("U_window_offset", 0)),
      fOffsetV(p.get<int>("V_window_offset", 0)),
      fOffsetX(p.get<int>("X_window_offset", 0)),
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
    mctruth_tree->Branch("Event", &fEventID);
    mctruth_tree->Branch("Run", &fRun);
    mctruth_tree->Branch("SubRun", &fSubRun);
    mctruth_tree->Branch("block_id", &mctruth_id);
    mctruth_tree->Branch("truth_track_id", &mctruth_trackid);
    mctruth_tree->Branch("pdg", &mctruth_pdg);
    mctruth_tree->Branch("generator_name", &mctruth_gen_name);
    mctruth_tree->Branch("status_code", &mctruth_status);
    mctruth_tree->Branch("x", &mctruth_x);
    mctruth_tree->Branch("y", &mctruth_y);
    mctruth_tree->Branch("z", &mctruth_z);
    mctruth_tree->Branch("px", &mctruth_Px);
    mctruth_tree->Branch("py", &mctruth_Py);
    mctruth_tree->Branch("pz", &mctruth_Pz);
    mctruth_tree->Branch("energy", &mctruth_en);
    mctruth_tree->Branch("process", &mctruth_process);

    mcneutrino_tree = tfs->make<TTree>("mcneutrinos", "mcneutrinos");
    mcneutrino_tree->Branch("Event", &fEventID);
    mcneutrino_tree->Branch("Run", &fRun);
    mcneutrino_tree->Branch("SubRun", &fSubRun);
    mcneutrino_tree->Branch("block_id", &mctruth_id);
    mcneutrino_tree->Branch("generator_name", &mctruth_gen_name);
    mcneutrino_tree->Branch("nupdg", &mcneutrino_nupdg);
    mcneutrino_tree->Branch("leptonpdg", &mcneutrino_leptonpdg);
    mcneutrino_tree->Branch("ccnc", &mcneutrino_ccnc);
    mcneutrino_tree->Branch("mode", &mcneutrino_mode);
    mcneutrino_tree->Branch("interactionType", &mcneutrino_iteractionType);
    mcneutrino_tree->Branch("target", &mcneutrino_target);
    mcneutrino_tree->Branch("hitnuc", &mcneutrino_hitnuc);
    mcneutrino_tree->Branch("hitquark", &mcneutrino_hitquark);
    mcneutrino_tree->Branch("w", &mcneutrino_w);
    mcneutrino_tree->Branch("x", &mcneutrino_x);
    mcneutrino_tree->Branch("y", &mcneutrino_y);
    mcneutrino_tree->Branch("qsqr", &mcneutrino_qsqr);
    mcneutrino_tree->Branch("pt", &mcneutrino_pt);
    mcneutrino_tree->Branch("theta", &mcneutrino_theta);
  }
  if (dump_mcparticles) {
    mcparticle_tree = tfs->make<TTree>("mcparticles", "mcparticles");
    mcparticle_tree->Branch("Event", &fEventID);
    mcparticle_tree->Branch("Run", &fRun);
    mcparticle_tree->Branch("SubRun", &fSubRun);
    mcparticle_tree->Branch("pdg", &mcparticle_pdg);
    mcparticle_tree->Branch("generator_name", &mcparticle_gen_name);
    mcparticle_tree->Branch("status_code", &mcparticle_status);
    mcparticle_tree->Branch("g4_track_id", &mcparticle_trackid);
    mcparticle_tree->Branch("mother", &mcparticle_mother);
    mcparticle_tree->Branch("truth_block_id", &mcparticle_truthid);
    mcparticle_tree->Branch("x", &mcparticle_x);
    mcparticle_tree->Branch("y", &mcparticle_y);
    mcparticle_tree->Branch("z", &mcparticle_z);
    mcparticle_tree->Branch("t", &mcparticle_t);
    mcparticle_tree->Branch("end_x", &mcparticle_end_x);
    mcparticle_tree->Branch("end_y", &mcparticle_end_y);
    mcparticle_tree->Branch("end_z", &mcparticle_end_z);
    mcparticle_tree->Branch("end_t", &mcparticle_end_t);
    mcparticle_tree->Branch("px", &mcparticle_Px);
    mcparticle_tree->Branch("py", &mcparticle_Py);
    mcparticle_tree->Branch("pz", &mcparticle_Pz);
    mcparticle_tree->Branch("energy", &mcparticle_en);
    mcparticle_tree->Branch("edep", &mcparticle_edep);
    mcparticle_tree->Branch("numelectrons", &mcparticle_numelectrons);
    mcparticle_tree->Branch("shower_edep", &mcparticle_shower_edep);
    mcparticle_tree->Branch("shower_numelectrons",
                            &mcparticle_shower_numelectrons);
    mcparticle_tree->Branch("process", &mcparticle_process);
  }
  if (dump_simides) {
    simide_tree = tfs->make<TTree>("simides", "simides");
    simide_tree->Branch("Event", &fEventID);
    simide_tree->Branch("Run", &fRun);
    simide_tree->Branch("SubRun", &fSubRun);
    simide_tree->Branch("channel", &sim_channel_id);
    simide_tree->Branch("timestamp", &tdc);
    simide_tree->Branch("numelectrons", &ide_numElectrons);
    simide_tree->Branch("energy", &ide_energy);
    simide_tree->Branch("x", &ide_x);
    simide_tree->Branch("y", &ide_y);
    simide_tree->Branch("z", &ide_z);
    simide_tree->Branch("trackID", &ide_trkId);
    simide_tree->Branch("origTrackID", &ide_origTrkId);
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
      mctruth_id = truth_block_counter;
      // NOTE: here we are making an assumption that the geant4 stage's process
      // name is largeant. This should be safe mostly.
      art::FindManyP<simb::MCParticle> assns(mctruthHandle, e, "largeant");
      for (size_t i = 0; i < mctruthHandle->size(); i++) {
        const simb::MCTruth &truthblock =
            *art::Ptr<simb::MCTruth>(mctruthHandle, i);
        if (truthblock.NeutrinoSet()) {
          const simb::MCNeutrino &mcneutrino = truthblock.GetNeutrino();
          mcneutrino_nupdg = mcneutrino.Nu().PdgCode();
          mcneutrino_leptonpdg = mcneutrino.Lepton().PdgCode();
          mcneutrino_ccnc = mcneutrino.CCNC();
          mcneutrino_mode = mcneutrino.Mode();
          mcneutrino_iteractionType = mcneutrino.InteractionType();
          mcneutrino_target = mcneutrino.Target();
          mcneutrino_hitnuc = mcneutrino.HitNuc();
          mcneutrino_hitquark = mcneutrino.HitQuark();
          mcneutrino_w = mcneutrino.W();
          mcneutrino_x = mcneutrino.X();
          mcneutrino_y = mcneutrino.Y();
          mcneutrino_qsqr = mcneutrino.QSqr();
          mcneutrino_pt = mcneutrino.Pt();
          mcneutrino_theta = mcneutrino.Theta();
          mcneutrino_tree->Fill();
        }
        int nparticles = truthblock.NParticles();
        if (dump_mcparticles) {
          std::vector<art::Ptr<simb::MCParticle>> matched_mcparts = assns.at(i);
          for (art::Ptr<simb::MCParticle> mcpart : matched_mcparts) {
            trkId_to_truthBlockId[mcpart->TrackId()] = truth_block_counter;
          }
        }
        for (int ipart = 0; ipart < nparticles; ipart++) {
          const simb::MCParticle &part = truthblock.GetParticle(ipart);
          mctruth_pdg = part.PdgCode();
          mctruth_gen_name = generator_name;
          mctruth_status = part.StatusCode();
          mctruth_process = part.Process();
          mctruth_trackid = part.TrackId();
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
  if (dump_simides) {
    auto simchannels =
        e.getValidHandle<std::vector<sim::SimChannel>>(simchannel_tag);
    track_en_sums.clear();
    track_electron_sums.clear();
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
          track_en_sums[ide.trackID] += ide.energy;
          track_electron_sums[ide.trackID] += ide.numElectrons;
          simide_tree->Fill();
        }
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
        mcparticle_mother = part.Mother();
        mcparticle_truthid =
            dump_mctruths ? trkId_to_truthBlockId.at(part.TrackId()) : -1;
        mcparticle_process = part.Process();
        mcparticle_x = part.Vx();
        mcparticle_y = part.Vy();
        mcparticle_z = part.Vz();
        mcparticle_t = part.T();
        mcparticle_end_x = part.EndX();
        mcparticle_end_y = part.EndY();
        mcparticle_end_z = part.EndZ();
        mcparticle_end_t = part.EndT();
        mcparticle_Px = part.Px();
        mcparticle_Py = part.Py();
        mcparticle_Pz = part.Pz();
        mcparticle_en = part.E();
        mcparticle_edep = track_en_sums.count(part.TrackId())
                              ? track_en_sums.at(part.TrackId())
                              : 0;
        mcparticle_numelectrons = track_electron_sums.count(part.TrackId())
                                      ? track_electron_sums.at(part.TrackId())
                                      : 0;
        mcparticle_shower_edep = track_en_sums.count(-part.TrackId())
                                     ? track_en_sums.at(-part.TrackId())
                                     : 0;
        mcparticle_shower_numelectrons =
            track_electron_sums.count(-part.TrackId())
                ? track_electron_sums.at(-part.TrackId())
                : 0;
        mcparticle_tree->Fill();
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
        tp_bufs[map_tag].from_tp(tp);
        ChannelInfo chinfo = get_channel_info_for_channel(geom, tp.channel);
        tp_channel_info_bufs[map_tag] = chinfo;
        if (tp_backtracking) {
          std::vector<sim::IDE> matched_ides =
              match_simides_to_tps(tp_bufs[map_tag], chinfo);
          tp_bufs[map_tag].populate_backtracking_info(matched_ides);
        }
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
      for (size_t i = 0; i < taHandle->size(); i++) {
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
            tp_bufs[map_tpInTaTag].from_tp(*tp);
            tp_channel_info_bufs[map_tag] =
                get_channel_info_for_channel(geom, tp->channel);
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
      for (size_t i = 0; i < tcHandle->size(); i++) {
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

void dunetrigger::TriggerAnaTree::make_tp_tree_if_needed(std::string tag,
                                                         bool assn) {
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
    TriggerPrimitiveBuffer &tp = tp_bufs[map_tag];
    tree->Branch("Event", &fEventID);
    tree->Branch("Run", &fRun);
    tree->Branch("SubRun", &fSubRun);
    tree->Branch("version", &tp.version);
    tree->Branch("flag", &tp.flag);
    tree->Branch("detid", &tp.detid);
    tree->Branch("channel", &tp.channel);
    tree->Branch("samples_over_threshold", &tp.samples_over_threshold);
    tree->Branch("time_start", &tp.time_start);
    tree->Branch("samples_to_peak", &tp.samples_to_peak);
    tree->Branch("adc_integral", &tp.adc_integral);
    tree->Branch("adc_peak", &tp.adc_peak);
    ChannelInfo &chinfo = tp_channel_info_bufs[map_tag];
    tree->Branch("ropid", &chinfo.rop_id);
    tree->Branch("view", &chinfo.view);
    tree->Branch("TPCSetID", &chinfo.tpcset_id);
    if (assn)
      tree->Branch("TAnumber", &fAssnIdx);
    if (tp_backtracking) {
      tree->Branch("bt_primary_track_id", &tp.bt_primary_track_id);
      tree->Branch("bt_primary_track_numelectron_frac",
                   &tp.bt_primary_track_numelectron_frac);
      tree->Branch("bt_primary_track_energy_frac",
                   &tp.bt_primary_track_energy_frac);
      tree->Branch("bt_edep", &tp.bt_edep);
      tree->Branch("bt_numelectrons", &tp.bt_numelectrons);
      tree->Branch("bt_x", &tp.bt_x);
      tree->Branch("bt_y", &tp.bt_y);
      tree->Branch("bt_z", &tp.bt_z);
      tree->Branch("bt_primary_x", &tp.bt_primary_x);
      tree->Branch("bt_primary_y", &tp.bt_primary_y);
      tree->Branch("bt_primary_z", &tp.bt_primary_z);
    }
  }
}

void dunetrigger::TriggerAnaTree::make_ta_tree_if_needed(std::string tag,
                                                         bool assn) {
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
    tree->Branch("Event", &fEventID);
    tree->Branch("Run", &fRun);
    tree->Branch("SubRun", &fSubRun);
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
    tree->Branch("type", &ta.type);
    tree->Branch("algorithm", &ta.algorithm);
    if (assn)
      tree->Branch("TCnumber", &fAssnIdx);
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
    tree->Branch("Event", &fEventID);
    tree->Branch("Run", &fRun);
    tree->Branch("SubRun", &fSubRun);
    tree->Branch("version", &tc.version);
    tree->Branch("time_start", &tc.time_start);
    tree->Branch("time_end", &tc.time_end);
    tree->Branch("time_candidate", &tc.time_candidate);
    tree->Branch("detid", &tc.detid);
    tree->Branch("type", &tc.type);
    tree->Branch("algorithm", &tc.algorithm);
  }
}

dunetrigger::ChannelInfo
dunetrigger::TriggerAnaTree::get_channel_info_for_channel(
    geo::WireReadoutGeom const *geom, int channel) {
  readout::ROPID rop = geom->ChannelToROP(channel);
  ChannelInfo result;
  result.rop_id = rop.ROP;
  result.tpcset_id = rop.asTPCsetID().TPCset;
  result.view = geom->View(rop);
  return result;
}

void dunetrigger::TriggerPrimitiveBuffer::from_tp(const TriggerPrimitive &tp) {
  version = tp.version;
  flag = 0;
  detid = tp.detid;
  channel = tp.channel;
  samples_over_threshold = (tp.time_over_threshold) /
                           dunetrigger::TPAlgTPCTool::ADC_SAMPLING_RATE_IN_DTS;
  time_start = tp.time_start;
  samples_to_peak = (tp.time_peak - tp.time_start) /
                    dunetrigger::TPAlgTPCTool::ADC_SAMPLING_RATE_IN_DTS;
  adc_integral = tp.adc_integral;
  adc_peak = tp.adc_peak;
}

std::vector<sim::IDE> dunetrigger::TriggerAnaTree::match_simides_to_tps(
    const TriggerPrimitiveBuffer &tp, const ChannelInfo &chinfo) const {
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  int offset = 0;
  switch (chinfo.view) {
  case geo::kU:
    offset += fOffsetU;
    break;
  case geo::kV:
    offset += fOffsetV;
    break;
  case geo::kW:
    offset += fOffsetX;
    break;
  default:
    break;
  }
  int sample_start = tp.time_start / TPAlgTPCTool::ADC_SAMPLING_RATE_IN_DTS;
  int sample_end = sample_start + tp.samples_over_threshold;
  sample_start += offset;
  sample_end += offset;
  sample_start = std::max(0, sample_start);
  sample_end = std::max(0, sample_end);

  if (sample_start > sample_end) {
    throw std::runtime_error("Invalid sample range");
  }
  art::Ptr<sim::SimChannel> sim_channel = bt_serv->FindSimChannel(tp.channel);
  std::vector<sim::IDE> matched_ides =
      sim_channel->TrackIDsAndEnergies(sample_start, sample_end);
  return matched_ides;
}

void dunetrigger::TriggerPrimitiveBuffer::populate_backtracking_info(
    const std::vector<sim::IDE> &ides) {
  bt_primary_track_id = INVALID;
  bt_primary_track_numelectron_frac = INVALID;
  bt_primary_track_energy_frac = INVALID;
  bt_edep = 0;
  bt_numelectrons = 0;
  bt_x = INVALID;
  bt_y = INVALID;
  bt_z = INVALID;
  bt_primary_x = INVALID;
  bt_primary_y = INVALID;
  bt_primary_z = INVALID;
  if (ides.size() == 0) {
    return;
  }
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  std::map<int, double> track_numelectrons;
  std::map<int, double> track_energies;
  for (const sim::IDE &ide : ides) {
    track_numelectrons[ide.trackID] += ide.numElectrons;
    bt_numelectrons += ide.numElectrons;
    bt_edep += ide.energy;
  }
  bt_primary_track_id =
      std::max_element(
          track_numelectrons.begin(), track_numelectrons.end(),
          [](const auto &a, const auto &b) { return a.second < b.second; })
          ->first;
  std::vector<sim::IDE> primary_ides;
  for (const sim::IDE &ide : ides) {
    if (ide.trackID == bt_primary_track_id) {
      primary_ides.push_back(ide);
    }
  }
  bt_primary_track_numelectron_frac =
      track_numelectrons[bt_primary_track_id] / bt_numelectrons;
  bt_primary_track_energy_frac = track_energies[bt_primary_track_id] / bt_edep;
  std::vector<double> bt_position = bt_serv->SimIDEsToXYZ(ides);
  std::vector<double> primary_bt_position = bt_serv->SimIDEsToXYZ(primary_ides);
  bt_x = bt_position[0];
  bt_y = bt_position[1];
  bt_z = bt_position[2];
  bt_primary_x = primary_bt_position[0];
  bt_primary_y = primary_bt_position[1];
  bt_primary_z = primary_bt_position[2];
}

DEFINE_ART_MODULE(dunetrigger::TriggerAnaTree)
