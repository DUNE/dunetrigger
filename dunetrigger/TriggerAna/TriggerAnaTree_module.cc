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

#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/WireReadout.h"
#include "larcoreobj/SimpleTypesAndConstants/readout_types.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include <regex>

#include <TDirectory.h>
#include <TFile.h>
#include <TTree.h>


#include <nlohmann/json.hpp>
using json = nlohmann::json;

using dunedaq::trgdataformats::TriggerActivityData;
using dunedaq::trgdataformats::TriggerCandidateData;
using dunedaq::trgdataformats::TriggerPrimitive;

// #define INVALID -99999
constexpr int INVALID=-99999;

namespace dunetrigger {
class TriggerAnaTree;

struct ChannelInfo {
  unsigned int rop_id;
  int view;
  unsigned int tpcset_id;
};

// Standard structure for all trees
struct EventDataBuffer {
  // buffer data members
  int event;
  int run;
  int subrun;

  void branch_on(TTree *tree) {
    tree->Branch("event", &event, "event/i");
    tree->Branch("run", &run, "run/i");
    tree->Branch("subrun", &subrun, "subrun/i");
  }

  void clear() {
    event = -1;
    run = -1;
    subrun = -1;
  }
};

struct EventSummaryBuffer {
  /**
   * Buffer of EventSummary data members
   */

  int mctruths_count;
  int mcparticles_count;
  int mcneutrinos_count;
  int simides_count;
  double tot_visible_energy_rop0, tot_visible_energy_rop1, tot_visible_energy_rop2, tot_visible_energy_rop3; // total visible energy per readout plane ID
  double tot_numelectrons_rop0, tot_numelectrons_rop1, tot_numelectrons_rop2, tot_numelectrons_rop3;

  void branch_on(TTree *tree) {
    // TODO: add num of tps for all collections?
    
    
    tree->Branch("mctruths_count", &mctruths_count);
    tree->Branch("mcparticles_count", &mcparticles_count);
    tree->Branch("mcneutrinos_count", &mcneutrinos_count);
    tree->Branch("simides_count", &simides_count);

    tree->Branch("tot_visible_energy_rop0", &tot_visible_energy_rop0);
    tree->Branch("tot_visible_energy_rop1", &tot_visible_energy_rop1);
    tree->Branch("tot_visible_energy_rop2", &tot_visible_energy_rop2);
    tree->Branch("tot_visible_energy_rop3", &tot_visible_energy_rop3);
    tree->Branch("tot_numelectrons_rop0", &tot_numelectrons_rop0);
    tree->Branch("tot_numelectrons_rop1", &tot_numelectrons_rop1);
    tree->Branch("tot_numelectrons_rop2", &tot_numelectrons_rop2);
    tree->Branch("tot_numelectrons_rop3", &tot_numelectrons_rop3);
  }


  void clear() {
    mctruths_count = -1;
    mcparticles_count = -1;
    mcneutrinos_count = -1;
    simides_count = -1;

    tot_visible_energy_rop0 = -1;
    tot_visible_energy_rop1 = -1;
    tot_visible_energy_rop2 = -1;
    tot_visible_energy_rop3 = -1;
  }
};

struct MCTruthBuffer {
  /**
   * Buffer of MCTruth data members
   */
  
  std::vector<int> pdg;
  std::vector<std::string> process;
  std::vector<int> status, id, trackid;
  std::vector<std::string> gen_name;
  std::vector<double> x, y, z, t;
  std::vector<double> Px, Py, Pz, P;
  std::vector<double> en, ek;

  void branch_on(TTree *tree) {
    tree->Branch("block_id", &id);
    tree->Branch("truth_track_id", &trackid);
    tree->Branch("pdg", &pdg);
    tree->Branch("generator_name", &gen_name);
    tree->Branch("status_code", &status);
    tree->Branch("x", &x);
    tree->Branch("y", &y);
    tree->Branch("z", &z);
    tree->Branch("t", &t);
    tree->Branch("px", &Px);
    tree->Branch("py", &Py);
    tree->Branch("pz", &Pz);
    tree->Branch("p", &P);
    tree->Branch("energy", &en);
    tree->Branch("kinetic_energy", &ek);
    tree->Branch("process", &process);
  }

  void clear() {
    pdg.clear();
    process.clear();
    status.clear();
    id.clear();
    trackid.clear();
    gen_name.clear();
    x.clear();
    y.clear();
    z.clear();
    t.clear();
    Px.clear();
    Py.clear();
    Pz.clear();
    P.clear();
    en.clear();
    ek.clear();
  }
};


struct MCNeutrinoBuffer {

  std::string mctruth_gen_name;
  int mctruth_id;
  int nupdg, leptonpdg, ccnc, mode, iteractionType,
      target, hitnuc, hitquark;
  double w, x, y, qsqr, pt, theta;

  void branch_on(TTree *tree) {
    tree->Branch("block_id", &mctruth_id);
    tree->Branch("generator_name", &mctruth_gen_name);
    tree->Branch("nupdg", &nupdg);
    tree->Branch("leptonpdg", &leptonpdg);
    tree->Branch("ccnc", &ccnc);
    tree->Branch("mode", &mode);
    tree->Branch("interactionType", &iteractionType);
    tree->Branch("target", &target);
    tree->Branch("hitnuc", &hitnuc);
    tree->Branch("hitquark", &hitquark);
    tree->Branch("w", &w);
    tree->Branch("x", &x);
    tree->Branch("y", &y);
    tree->Branch("qsqr", &qsqr);
    tree->Branch("pt", &pt);
    tree->Branch("theta", &theta);
  }
};


struct MCParticleBuffer {
  /**
   * Buffer of MCParticles data members
   */
  int pdg;
  std::string process;
  int status, trackid, truthid, mother;
  std::string gen_name;
  double x, y, z, t;
  double end_x, end_y, end_z, end_t;
  double Px, Py, Pz;
  double en, ek;
  double edep, numelectrons;
  double shower_edep, shower_numelectrons;

  void branch_on(TTree *tree) {
    tree->Branch("pdg", &pdg);
    tree->Branch("generator_name", &gen_name);
    tree->Branch("status_code", &status);
    tree->Branch("g4_track_id", &trackid);
    tree->Branch("mother", &mother);
    tree->Branch("truth_block_id", &truthid);
    tree->Branch("x", &x);
    tree->Branch("y", &y);
    tree->Branch("z", &z);
    tree->Branch("t", &t);
    tree->Branch("end_x", &end_x);
    tree->Branch("end_y", &end_y);
    tree->Branch("end_z", &end_z);
    tree->Branch("end_t", &end_t);
    tree->Branch("px", &Px);
    tree->Branch("py", &Py);
    tree->Branch("pz", &Pz);
    tree->Branch("energy", &en);
    tree->Branch("kinetic_energy", &ek);
    tree->Branch("edep", &edep);
    tree->Branch("numelectrons", &numelectrons);
    tree->Branch("shower_edep", &shower_edep);
    tree->Branch("shower_numelectrons", &shower_numelectrons);
    tree->Branch("process", &process);
  }

  void clear() {
    // FIXME: implement
  }
};


struct SimIDEBuffer {
  /**
   * Buffer of SimIDEBuffer data members
   */

  unsigned int sim_channel_id;
  int tdc;
  float numElectrons;
  float energy;
  float x;
  float y;
  float z;
  int trkId;
  float origTrkId;
  float readout_plane_id;
  float readout_view;
  float detector_element;
  // int particle_pdg;
  // int parent_pdg;

  void branch_on(TTree *tree) {
    tree->Branch("channel", &sim_channel_id);
    tree->Branch("timestamp", &tdc);
    tree->Branch("numelectrons", &numElectrons);
    tree->Branch("energy", &energy);
    tree->Branch("x", &x);
    tree->Branch("y", &y);
    tree->Branch("z", &z);
    tree->Branch("trackID", &trkId);
    tree->Branch("origTrackID", &origTrkId);

    tree->Branch("readout_plane_id", &readout_plane_id);
    tree->Branch("readout_view", &readout_view);
    tree->Branch("detector_element", &detector_element);

    // tree->Branch("pdg_id", &particle_pdg);
    // tree->Branch("parent_pdg_id", &parent_pdg);
  }

  void clear() {
    sim_channel_id = INVALID;
    tdc = INVALID;
    numElectrons = INVALID;
    energy = INVALID;
    x = INVALID;
    y = INVALID;
    z = INVALID;
    trkId = INVALID;
    origTrkId = INVALID;
    readout_plane_id = INVALID;
    readout_view = INVALID;
    detector_element = INVALID;
    // particle_pdg = INVALID;
    // parent_pdg = INVALID;
  }

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
  int bt_num_tracks;
  double bt_primary_track_numelectron_frac;
  double bt_primary_track_energy_frac;
  double bt_edep;
  double bt_numelectrons;
  double bt_x, bt_y, bt_z;
  double bt_primary_x, bt_primary_y, bt_primary_z;
  int bt_mctruth_block_id;
  std::string bt_mctruth_gen_name;

  // Adding chinfo here -- no need to track it separately
  ChannelInfo chinfo;

  // Populate the buffer from a TP object
  void from_tp(const TriggerPrimitive &tp);

  // Populate backtracking information matched ides
  void populate_backtracking_info(const std::vector<sim::IDE> &ides,
                                  const std::unordered_map<int, int> &trkid_to_truth_block,
                                  std::unordered_map<int, std::string> &truth_id_to_gen);

  void branch_on(TTree *tree, bool backtracking) {
    tree->Branch("version", &this->version);
    tree->Branch("flag", &this->flag);
    tree->Branch("detid", &this->detid);
    tree->Branch("channel", &this->channel);
    tree->Branch("samples_over_threshold", &this->samples_over_threshold);
    tree->Branch("time_start", &this->time_start);
    tree->Branch("samples_to_peak", &this->samples_to_peak);
    tree->Branch("adc_integral", &this->adc_integral);
    tree->Branch("adc_peak", &this->adc_peak);

    tree->Branch("readout_plane_id", &this->chinfo.rop_id);
    tree->Branch("readout_view", &this->chinfo.view);
    tree->Branch("TPCSetID", &this->chinfo.tpcset_id);

    // Add backtracking only if requested
    if (backtracking) {
      tree->Branch("bt_primary_track_id", &this->bt_primary_track_id);
      tree->Branch("bt_primary_track_numelectron_frac", &this->bt_primary_track_numelectron_frac);
      tree->Branch("bt_primary_track_energy_frac", &this->bt_primary_track_energy_frac);
      tree->Branch("bt_edep", &this->bt_edep);
      tree->Branch("bt_numelectrons", &this->bt_numelectrons);
      tree->Branch("bt_x", &this->bt_x);
      tree->Branch("bt_y", &this->bt_y);
      tree->Branch("bt_z", &this->bt_z);
      tree->Branch("bt_primary_x", &this->bt_primary_x);
      tree->Branch("bt_primary_y", &this->bt_primary_y);
      tree->Branch("bt_primary_z", &this->bt_primary_z);
      tree->Branch("bt_truth_block_id", &this->bt_mctruth_block_id);
      tree->Branch("bt_generator_name", &this->bt_mctruth_gen_name);
    }
  }

  // Reset the buffer
  void reset() {

    version = 0;
    flag = 0x0;
    detid = 0;

    channel = 0;
    samples_over_threshold = 0;
    time_start = 0;
    samples_to_peak = 0;
    adc_integral = 0;
    adc_peak = 0;

    // Reset backtracker members
    bt_primary_track_id = INVALID;
    bt_num_tracks = INVALID;
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
    bt_mctruth_block_id = INVALID;
    bt_mctruth_gen_name.clear();
  };
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
  void endJob() override;

private:
  art::ServiceHandle<art::TFileService> tfs;
  std::map<std::string, TTree *> tree_map;
  // buffers for writing to ROOT Trees
  EventDataBuffer ev_buf;

  size_t fAssnIdx;

  std::unordered_map<int, int> trkId_to_truthBlockId;
  std::unordered_map<int, std::string> truthBlockId_to_generator_name;
  std::map<std::string, TriggerPrimitiveBuffer> tp_bufs;
  std::map<std::string, ChannelInfo> tp_channel_info_bufs;
  std::map<std::string, TriggerActivityData> ta_bufs;
  std::map<std::string, TriggerCandidateData> tc_bufs;
  std::map<int, double> track_en_sums;
  std::map<int, double> track_electron_sums;

  bool dump_tp, dump_ta, dump_tc;
  std::string tp_tag_regex, ta_tag_regex, tc_tag_regex;
  

  bool tp_backtracking;

  void make_tp_tree_if_needed(std::string tag, bool assn = false);
  void make_ta_tree_if_needed(std::string tag, bool assn = false);
  void make_tc_tree_if_needed(std::string tag);

  std::vector<sim::IDE> match_simides_to_tps(const TriggerPrimitiveBuffer &tp, const std::string &tool_type) const;

  ChannelInfo get_channel_info_for_channel(geo::WireReadoutGeom const *geom, int channel);

  bool dump_mctruths;
  TTree *mctruth_tree;
  MCTruthBuffer mctruth_buf;

  TTree *mcneutrino_tree;
  MCNeutrinoBuffer mcneutrino_buf;

  bool dump_mcparticles;
  TTree *mcparticle_tree;
  MCParticleBuffer mcparticle_buf;

  bool dump_summary_info;
  // visible energy for the event
  TTree *summary_tree;
  EventSummaryBuffer evsum_buf;

  std::map<std::string, std::array<int, 3>> bt_view_offsets;

  bool dump_simides;
  std::string simchannel_tag;
  TTree *simide_tree;
  SimIDEBuffer simide_buf;
  
  // JSON metadata
  json info_data;
  bool first_event_flag;
};

dunetrigger::TriggerAnaTree::TriggerAnaTree(fhicl::ParameterSet const &p)
    : EDAnalyzer{p}, 
    dump_tp(p.get<bool>("dump_tp")),
    dump_ta(p.get<bool>("dump_ta")),
    dump_tc(p.get<bool>("dump_tc")),
    tp_tag_regex(p.get<std::string>("tp_tag_regex", ".*")),
    ta_tag_regex(p.get<std::string>("ta_tag_regex", ".*")),
    tc_tag_regex(p.get<std::string>("tc_tag_regex", ".*")),
    tp_backtracking(p.get<bool>("tp_backtracking", false)),
    dump_mctruths(p.get<bool>("dump_mctruths", true)),
    dump_mcparticles(p.get<bool>("dump_mcparticles", true)),
    dump_summary_info(p.get<bool>("dump_summary_info", true)),
    dump_simides(p.get<bool>("dump_simides", true)),
    simchannel_tag(p.get<std::string>("simchannel_tag", "tpcrawdecoder:simpleSC"))
// More initializers here.
{
  std::vector<fhicl::ParameterSet> offsets = p.get<std::vector<fhicl::ParameterSet>>("bt_window_offsets");
  for (const auto &offset : offsets) {
    bt_view_offsets[offset.get<std::string>("tool_type")] = {offset.get<int>("U"), offset.get<int>("V"),
                                                             offset.get<int>("X")};
  }
}

void dunetrigger::TriggerAnaTree::beginJob() {
  if (dump_mctruths) {
    mctruth_tree = tfs->make<TTree>("mctruths", "mctruths");
    ev_buf.branch_on(mctruth_tree);
    mctruth_buf.branch_on(mctruth_tree);

    mcneutrino_tree = tfs->make<TTree>("mcneutrinos", "mcneutrinos");
    ev_buf.branch_on(mcneutrino_tree);
    mcparticle_buf.branch_on(mcneutrino_tree);

  }

  if (dump_mcparticles) {
    mcparticle_tree = tfs->make<TTree>("mcparticles", "mcparticles");
    ev_buf.branch_on(mcparticle_tree);
    mcparticle_buf.branch_on(mcparticle_tree);
  }

  if (dump_simides) {
    simide_tree = tfs->make<TTree>("simides", "simides");
    ev_buf.branch_on(simide_tree);
    simide_buf.branch_on(simide_tree);
  }

  if (dump_summary_info) {
    summary_tree = tfs->make<TTree>("event_summary", "event_summary");
    ev_buf.branch_on(summary_tree);
    evsum_buf.branch_on(summary_tree);
  }

  // Save detector settings
  auto const &geo = art::ServiceHandle<geo::Geometry>();

  // Geometry
  info_data["geo"] = {};
  info_data["geo"]["detector"] = geo->DetectorName();

  if (tp_backtracking) {
    for (const auto &[tool, offsets] : bt_view_offsets) {
      info_data["backtracker"][tool]["offset_U"] = offsets[0];
      info_data["backtracker"][tool]["offset_V"] = offsets[1];
      info_data["backtracker"][tool]["offset_X"] = offsets[2];
    }
  }

  first_event_flag = true;
}

void dunetrigger::TriggerAnaTree::analyze(art::Event const &e) {

  ev_buf.run = e.run();
  ev_buf.subrun = e.subRun();
  ev_buf.event = e.event();

  // reset visible energy counters
  mctruth_buf.clear();
  evsum_buf.clear();

  size_t mctruths_count = 0;
  size_t mcneutrinos_count = 0;
  size_t mcparticles_count = 0;
  size_t simides_count = 0;

  // get a service handle for geometry
  geo::WireReadoutGeom const *geom = &art::ServiceHandle<geo::WireReadout>()->Get();

  if (dump_mctruths) {
    std::vector<art::Handle<std::vector<simb::MCTruth>>> mctruthHandles = e.getMany<std::vector<simb::MCTruth>>();

    int truth_block_counter = 0;
    trkId_to_truthBlockId.clear();
    truthBlockId_to_generator_name.clear();

    for (auto const &mctruthHandle : mctruthHandles) {
      // Extract the generator name from the truth handle input label
      std::string generator_name = mctruthHandle.provenance()->inputTag().label();
      // Store generator name for TP backtracking
      truthBlockId_to_generator_name[truth_block_counter] = generator_name;

      mctruth_buf.id.push_back(truth_block_counter);
      // NOTE: here we are making an assumption that the geant4 stage's process
      // name is largeant. This should be safe mostly.
      art::FindManyP<simb::MCParticle> assns(mctruthHandle, e, "largeant");
      for (size_t i = 0; i < mctruthHandle->size(); i++) {
        const simb::MCTruth &truthblock = *art::Ptr<simb::MCTruth>(mctruthHandle, i);


        std::vector<art::Ptr<simb::MCParticle>> matched_mcparts = assns.at(i);
        for (art::Ptr<simb::MCParticle> mcpart : matched_mcparts) {
          trkId_to_truthBlockId[mcpart->TrackId()] = truth_block_counter;
        }

        if (truthblock.NeutrinoSet()) {

          const simb::MCNeutrino &mcneutrino = truthblock.GetNeutrino();
          
          mcneutrino_buf.mctruth_id = truth_block_counter;
          mcneutrino_buf.mctruth_gen_name = generator_name;
          mcneutrino_buf.nupdg = mcneutrino.Nu().PdgCode();
          mcneutrino_buf.leptonpdg = mcneutrino.Lepton().PdgCode();
          mcneutrino_buf.ccnc = mcneutrino.CCNC();
          mcneutrino_buf.mode = mcneutrino.Mode();
          mcneutrino_buf.iteractionType = mcneutrino.InteractionType();
          mcneutrino_buf.target = mcneutrino.Target();
          mcneutrino_buf.hitnuc = mcneutrino.HitNuc();
          mcneutrino_buf.hitquark = mcneutrino.HitQuark();
          mcneutrino_buf.w = mcneutrino.W();
          mcneutrino_buf.x = mcneutrino.X();
          mcneutrino_buf.y = mcneutrino.Y();
          mcneutrino_buf.qsqr = mcneutrino.QSqr();
          mcneutrino_buf.pt = mcneutrino.Pt();
          mcneutrino_buf.theta = mcneutrino.Theta();
          mcneutrino_tree->Fill();
          ++mcneutrinos_count;
        }

        int nparticles = truthblock.NParticles();


        for (int ipart = 0; ipart < nparticles; ipart++) {

          const simb::MCParticle &part = truthblock.GetParticle(ipart);
          mctruth_buf.id.push_back(truth_block_counter);
          mctruth_buf.pdg.push_back(part.PdgCode());
          mctruth_buf.gen_name.push_back(generator_name);
          mctruth_buf.status.push_back(part.StatusCode());
          mctruth_buf.process.push_back(part.Process());
          mctruth_buf.trackid.push_back(part.TrackId());
          mctruth_buf.x.push_back(part.Vx());
          mctruth_buf.y.push_back(part.Vy());
          mctruth_buf.z.push_back(part.Vz());
          mctruth_buf.t.push_back(part.T());
          mctruth_buf.Px.push_back(part.Px());
          mctruth_buf.Py.push_back(part.Py());
          mctruth_buf.Pz.push_back(part.Pz());
          mctruth_buf.P.push_back(part.P());
          mctruth_buf.en.push_back(part.E());
          mctruth_buf.ek.push_back(part.E() - part.Mass());
          mctruth_tree->Fill();
          ++mctruths_count;
        }
        truth_block_counter++;
      }
    }

    json j_mctruth_gen_map(truthBlockId_to_generator_name);
    info_data["mcthruth_blockid_map"] = j_mctruth_gen_map;

  }

  if (dump_simides) {
    auto simchannels = e.getValidHandle<std::vector<sim::SimChannel>>(simchannel_tag);
    art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;

    track_en_sums.clear();
    track_electron_sums.clear();
    for (const sim::SimChannel &sc : *simchannels) {
      simide_buf.sim_channel_id = sc.Channel();

      sim::SimChannel::TDCIDEs_t const &tdcidemap = sc.TDCIDEMap();
      for (const sim::TDCIDE &tdcide : tdcidemap) {
        simide_buf.tdc = tdcide.first;
        for (const sim::IDE& ide : tdcide.second) {

          simide_buf.numElectrons = ide.numElectrons;
          simide_buf.energy = ide.energy;
          simide_buf.x = ide.x;
          simide_buf.y = ide.y;
          simide_buf.z = ide.z;
          simide_buf.trkId = ide.trackID;
          simide_buf.origTrkId = ide.origTrackID;


          track_en_sums[ide.trackID] += ide.energy;
          track_electron_sums[ide.trackID] += ide.numElectrons;
          // higher level geometric info for IDEs
          ChannelInfo chinfo = get_channel_info_for_channel(geom, simide_buf.sim_channel_id);
          simide_buf.readout_plane_id = chinfo.rop_id;
          simide_buf.readout_view = chinfo.view;
          simide_buf.detector_element = chinfo.tpcset_id; // APA/CRP ID
          
          // populate the total visible energy counters by plane
          if (chinfo.rop_id == 0) {
            evsum_buf.tot_visible_energy_rop0 += ide.energy;
            evsum_buf.tot_numelectrons_rop0 += ide.numElectrons;
          }
          else if (chinfo.rop_id == 1) {
            evsum_buf.tot_visible_energy_rop1 += ide.energy;
            evsum_buf.tot_numelectrons_rop1 += ide.numElectrons;
          }
          else if (chinfo.rop_id == 2) {
            evsum_buf.tot_visible_energy_rop2 += ide.energy;
            evsum_buf.tot_numelectrons_rop2 += ide.numElectrons;
          }
          else if (chinfo.rop_id == 3) {
            evsum_buf.tot_visible_energy_rop3 += ide.energy;
            evsum_buf.tot_numelectrons_rop3 += ide.numElectrons;
          }
          simide_tree->Fill();
          ++simides_count;
        }
      }
    }
  }

  if (dump_mcparticles) {

    std::vector<art::Handle<std::vector<simb::MCParticle>>> mcparticleHandles =
        e.getMany<std::vector<simb::MCParticle>>();

    for (auto const &mcparticleHandle : mcparticleHandles) {

      std::string generator_name = mcparticleHandle.provenance()->inputTag().label();

      for (const simb::MCParticle &part : *mcparticleHandle) {
        mcparticle_buf.pdg = part.PdgCode();
        mcparticle_buf.gen_name = generator_name;
        mcparticle_buf.status = part.StatusCode();
        mcparticle_buf.trackid = part.TrackId();
        mcparticle_buf.mother = part.Mother();
        mcparticle_buf.truthid = dump_mctruths ? trkId_to_truthBlockId.at(part.TrackId()) : -1;
        mcparticle_buf.process = part.Process();
        mcparticle_buf.x = part.Vx();
        mcparticle_buf.y = part.Vy();
        mcparticle_buf.z = part.Vz();
        mcparticle_buf.t = part.T();
        mcparticle_buf.end_x = part.EndX();
        mcparticle_buf.end_y = part.EndY();
        mcparticle_buf.end_z = part.EndZ();
        mcparticle_buf.end_t = part.EndT();
        mcparticle_buf.Px = part.Px();
        mcparticle_buf.Py = part.Py();
        mcparticle_buf.Pz = part.Pz();
        mcparticle_buf.en = part.E();
        mcparticle_buf.ek = part.E() - part.Mass();
        mcparticle_buf.edep = track_en_sums.count(part.TrackId()) ? track_en_sums.at(part.TrackId()) : 0;
        mcparticle_buf.numelectrons =
            track_electron_sums.count(part.TrackId()) ? track_electron_sums.at(part.TrackId()) : 0;
        mcparticle_buf.shower_edep = track_en_sums.count(-part.TrackId()) ? track_en_sums.at(-part.TrackId()) : 0;
        mcparticle_buf.shower_numelectrons =
            track_electron_sums.count(-part.TrackId()) ? track_electron_sums.at(-part.TrackId()) : 0;
        mcparticle_tree->Fill();
        ++mcparticles_count;
      }
    }
  }

  if (dump_tp) {
    std::vector<art::Handle<std::vector<TriggerPrimitive>>> tpHandles = e.getMany<std::vector<TriggerPrimitive>>();

    if ( first_event_flag ) {
      info_data["tpg"] = {};
    }

    std::regex tp_regex(this->tp_tag_regex);
    for (auto const &tpHandle : tpHandles) {

      std::string tag = tpHandle.provenance()->inputTag().encode();
      if ( !std::regex_match(tag, tp_regex) ) {
        continue;
      }

      fhicl::ParameterSet tp_params = tpHandle.provenance()->parameterSet().get<fhicl::ParameterSet>("tpalg");
      std::string tp_tool_type = tp_params.get<std::string>("tool_type");

      if ( first_event_flag ) {
        info_data["tpg"][tag]["tool"] = tp_tool_type;
        info_data["tpg"][tag]["threshold_tpg_plane0"] = tp_params.get<int>("threshold_tpg_plane0");
        info_data["tpg"][tag]["threshold_tpg_plane1"] = tp_params.get<int>("threshold_tpg_plane1");
        info_data["tpg"][tag]["threshold_tpg_plane2"] = tp_params.get<int>("threshold_tpg_plane2");
      }

      std::string map_tag = "tp/" + tag;

      TriggerPrimitiveBuffer &curr_tp_buf = tp_bufs[map_tag];
      make_tp_tree_if_needed(tag);

      TTree *cur_tp_tree = tree_map[map_tag];
      for (const TriggerPrimitive &tp : *tpHandle) {
        curr_tp_buf.from_tp(tp);
        curr_tp_buf.chinfo = get_channel_info_for_channel(geom, tp.channel);
        if (tp_backtracking) {
          std::vector<sim::IDE> matched_ides = match_simides_to_tps(curr_tp_buf, tp_tool_type);
          curr_tp_buf.populate_backtracking_info(matched_ides, trkId_to_truthBlockId, truthBlockId_to_generator_name);
        }
        cur_tp_tree->Fill();
      }
    }


    evsum_buf.mctruths_count = mctruths_count;
    evsum_buf.mcparticles_count = mcparticles_count;
    evsum_buf.mcneutrinos_count = mcneutrinos_count;
    evsum_buf.simides_count = simides_count;

  }

  if (dump_ta) {
    std::vector<art::Handle<std::vector<TriggerActivityData>>> taHandles =
        e.getMany<std::vector<TriggerActivityData>>();

    std::regex ta_regex(this->ta_tag_regex);

    for (auto const &taHandle : taHandles) {

      art::FindManyP<TriggerPrimitive> assns(taHandle, e, taHandle.provenance()->moduleLabel());
      std::string tag = taHandle.provenance()->inputTag().encode();
      if ( !std::regex_match(tag, ta_regex) ) {
        continue;
      }
      std::string map_tag = "ta/" + tag;
      make_ta_tree_if_needed(tag);
      for (size_t i = 0; i < taHandle->size(); i++) {
        const TriggerActivityData &ta = *art::Ptr<TriggerActivityData>(taHandle, i);
        if (assns.isValid()) {
          art::InputTag ta_input_tag = taHandle.provenance()->inputTag();
          std::string tpInTaTag =
              art::InputTag(ta_input_tag.label(), ta_input_tag.instance() + "inTAs", ta_input_tag.process()).encode();
          std::string map_tpInTaTag = "tp/" + tpInTaTag;
          make_tp_tree_if_needed(tpInTaTag, true);
          fAssnIdx = i;
          std::vector<art::Ptr<TriggerPrimitive>> matched_tps = assns.at(i);
          TriggerPrimitiveBuffer &curr_tp_buf = tp_bufs[map_tpInTaTag];
          TTree *cur_tp_tree = tree_map[map_tpInTaTag];

          for (art::Ptr<TriggerPrimitive> tp : matched_tps) {
            curr_tp_buf.from_tp(*tp);
            curr_tp_buf.chinfo = get_channel_info_for_channel(geom, tp->channel);
            cur_tp_tree->Fill();
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

    std::regex tc_regex(this->tc_tag_regex);

    for (auto const &tcHandle : tcHandles) {
      art::FindManyP<TriggerActivityData> assns(tcHandle, e, tcHandle.provenance()->moduleLabel());
      std::string tag = tcHandle.provenance()->inputTag().encode();
      if ( !std::regex_match(tag, tc_regex) ) {
        continue;
      }
      std::string map_tag = "tc/" + tag;
      make_tc_tree_if_needed(tag);
      for (size_t i = 0; i < tcHandle->size(); i++) {
        const TriggerCandidateData &tc = *art::Ptr<TriggerCandidateData>(tcHandle, i);
        if (assns.isValid()) {
          art::InputTag tc_input_tag = tcHandle.provenance()->inputTag();
          std::string taInTcTag =
              art::InputTag(tc_input_tag.label(), tc_input_tag.instance() + "inTCs", tc_input_tag.process()).encode();
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

  if (dump_summary_info) {
    summary_tree->Fill();
  }

  first_event_flag = false;
}

void dunetrigger::TriggerAnaTree::endJob() {

  art::ServiceHandle<art::TFileService> tfs;
  auto n = tfs->make<TNamed>("info", info_data.dump().c_str());
  n->Write();
}

void dunetrigger::TriggerAnaTree::make_tp_tree_if_needed(std::string tag, bool assn) {
  std::string map_tag = "tp/" + tag;
  if (!tree_map.count(map_tag)) {
    art::TFileDirectory tp_dir = tfs->mkdir("TriggerPrimitives", "Trigger Primitive Trees");
    std::cout << "Creating new TTree for " << tag << std::endl;

    // Replace ":" with "_" in TTree names so that they can be used in ROOT's
    // intepreter
    std::string tree_name = tag;
    std::replace(tree_name.begin(), tree_name.end(), ':', '_');

    // Create tree
    TTree *tree = tp_dir.make<TTree>(tree_name.c_str(), tree_name.c_str());
    tree_map[map_tag] = tree;

    // Initialize TP buffer
    TriggerPrimitiveBuffer &tp = tp_bufs[map_tag];

    ev_buf.branch_on(tree);
    tp.branch_on(tree, tp_backtracking);
    if (assn)
      tree->Branch("TAnumber", &fAssnIdx);
  }
}

void dunetrigger::TriggerAnaTree::make_ta_tree_if_needed(std::string tag, bool assn) {
  std::string map_tag = "ta/" + tag;
  if (!tree_map.count(map_tag)) {
    art::TFileDirectory ta_dir = tfs->mkdir("TriggerActivities", "Trigger Activity Trees");
    std::cout << "Creating new TTree for " << tag << std::endl;
    // Replace ":" with "_" in TTree names so that they can be used in ROOT's
    // intepreter
    std::string tree_name = tag;
    std::replace(tree_name.begin(), tree_name.end(), ':', '_');
    TTree *tree = ta_dir.make<TTree>(tree_name.c_str(), tree_name.c_str());
    tree_map[map_tag] = tree;
    TriggerActivityData &ta = ta_bufs[map_tag];
    ev_buf.branch_on(tree);

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
    art::TFileDirectory tc_dir = tfs->mkdir("TriggerCandidates", "Trigger Candidate Trees");
    std::cout << "Creating new TTree for " << tag << std::endl;
    // Replace ":" with "_" in TTree names so that they can be used in ROOT's
    // intepreter
    std::string tree_name = tag;
    std::replace(tree_name.begin(), tree_name.end(), ':', '_');
    TTree *tree = tc_dir.make<TTree>(tree_name.c_str(), tree_name.c_str());
    tree_map[map_tag] = tree;
    TriggerCandidateData &tc = tc_bufs[map_tag];
    ev_buf.branch_on(tree);

    tree->Branch("version", &tc.version);
    tree->Branch("time_start", &tc.time_start);
    tree->Branch("time_end", &tc.time_end);
    tree->Branch("time_candidate", &tc.time_candidate);
    tree->Branch("detid", &tc.detid);
    tree->Branch("type", &tc.type);
    tree->Branch("algorithm", &tc.algorithm);
  }
}

dunetrigger::ChannelInfo dunetrigger::TriggerAnaTree::get_channel_info_for_channel(geo::WireReadoutGeom const *geom,
                                                                                   int channel) {
  readout::ROPID rop = geom->ChannelToROP(channel);
  ChannelInfo result;
  result.rop_id = rop.ROP;
  result.tpcset_id = rop.asTPCsetID().TPCset;
  result.view = geom->View(rop);
  return result;
}

void dunetrigger::TriggerPrimitiveBuffer::from_tp(const TriggerPrimitive &tp) {
  version = 2; // temp, since variables below are converted to v2 version while the TP version in TriggerSim is still 1.
               // Go back to "= tp.version" after changing triggeralgs to v5 (and using TriggerPrimitive2.hpp as header)
  flag = 0;
  detid = tp.detid;
  channel = tp.channel;
  samples_over_threshold = (tp.time_over_threshold) / dunetrigger::TPAlgTPCTool::ADC_SAMPLING_RATE_IN_DTS;
  time_start = tp.time_start;
  samples_to_peak = (tp.time_peak - tp.time_start) / dunetrigger::TPAlgTPCTool::ADC_SAMPLING_RATE_IN_DTS;
  adc_integral = tp.adc_integral;
  adc_peak = tp.adc_peak;
}

void dunetrigger::TriggerPrimitiveBuffer::populate_backtracking_info(
    const std::vector<sim::IDE> &ides, 
    const std::unordered_map<int, int> &trkid_to_truth_block,
    std::unordered_map<int, std::string> &truth_id_to_gen) {

  // Reset backtracker members
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
  bt_mctruth_block_id = INVALID;
  bt_mctruth_gen_name.clear();

  if (ides.size() == 0) {
    return;
  }
  
  // Backtracked ides - consider only ides with valid track ID.
  std::vector<sim::IDE> bt_ides;
  bt_ides.reserve(ides.capacity());
  std::copy_if (ides.begin(), ides.end(), std::back_inserter(bt_ides), [](const sim::IDE& ide){return ide.trackID != 0;} );

  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;

  std::map<int, double> track_numelectrons;
  std::map<int, double> track_energies;

  if (bt_ides.empty() ) {
    std::cout << "Empty IDEs set!" << std::endl;
    return;
  }

  for (const sim::IDE &ide : bt_ides) {
    int mc_track_id = pi_serv->TrackIdToParticle_P(ide.trackID)->TrackId();
    track_numelectrons[mc_track_id] += ide.numElectrons;
    bt_numelectrons += ide.numElectrons;
    bt_edep += ide.energy;
  }
  bt_primary_track_id =
      std::max_element(track_numelectrons.begin(), track_numelectrons.end(), [](const auto &a, const auto &b) {
        return a.second < b.second;
      })->first;

  std::vector<sim::IDE> primary_ides;
  for (const sim::IDE &ide : bt_ides) {
    if (pi_serv->TrackIdToParticle_P(ide.trackID)->TrackId() == bt_primary_track_id) {
      primary_ides.push_back(ide);
    }
  }

  bt_num_tracks = track_numelectrons.size();
  bt_primary_track_numelectron_frac = track_numelectrons[bt_primary_track_id] / bt_numelectrons;
  bt_primary_track_energy_frac = track_energies[bt_primary_track_id] / bt_edep;

  std::vector<double> bt_position = bt_serv->SimIDEsToXYZ(ides);
  std::vector<double> primary_bt_position = bt_serv->SimIDEsToXYZ(primary_ides);

  bt_x = bt_position[0];
  bt_y = bt_position[1];
  bt_z = bt_position[2];
  bt_primary_x = primary_bt_position[0];
  bt_primary_y = primary_bt_position[1];
  bt_primary_z = primary_bt_position[2];

  bt_mctruth_block_id = trkid_to_truth_block.at(bt_primary_track_id);
  bt_mctruth_gen_name = truth_id_to_gen.at(bt_mctruth_block_id);
}


std::vector<sim::IDE> dunetrigger::TriggerAnaTree::match_simides_to_tps(const TriggerPrimitiveBuffer &tp,
                                                                        const std::string &tool_type) const {

  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  auto it = bt_view_offsets.find(tool_type);
  if (it == bt_view_offsets.end()) {
    std::cout << "No offsets found for tool type " << tool_type << ", using 0,0,0" << std::endl;
  }
  const std::array<int, 3> &offsets = it != bt_view_offsets.end() ? it->second : std::array<int, 3>{0, 0, 0};
  int offset = 0;
  switch (tp.chinfo.view) {
    case geo::kU:
      offset = offsets[0];
      break;
    case geo::kV:
      offset = offsets[1];
      break;
    case geo::kW:
      offset = offsets[2];
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
  std::vector<sim::IDE> matched_ides = sim_channel->TrackIDsAndEnergies(sample_start, sample_end);
  return matched_ides;
}

DEFINE_ART_MODULE(dunetrigger::TriggerAnaTree)
