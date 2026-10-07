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
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"

#include "detdataformats/trigger/TriggerActivityData.hpp"
#include "detdataformats/trigger/TriggerCandidateData.hpp"

#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/WireReadout.h"
#include "larcoreobj/SimpleTypesAndConstants/readout_types.h"
#include "lardataobj/Simulation/SimChannel.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include <regex>

#include "messagefacility/MessageLogger/MessageLogger.h"

#include <TNamed.h>
#include <TTree.h>


#include "VectorFieldsBuffer.hh"
#include "ScalarFieldsBuffer.hh"
#include "TriggerAnaTree_module.hh"

#include "dunetrigger/TriggerSim/TPAlgTools/TPAlgTPCTool.hh"
#include "dunetrigger/vendor/lardata/ArtDataHelper/GetManyByRegexTag.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include <algorithm>
#include <iostream>
#include <map>
#include <set>

#include <boost/range/adaptor/map.hpp>
#include <boost/range/algorithm/set_algorithm.hpp>


#include <nlohmann/json.hpp>
using json = nlohmann::json;

using dunedaq::trgdataformats::TriggerPrimitive;
using dunedaq::trgdataformats::TriggerActivityData;
using dunedaq::trgdataformats::TriggerCandidateData;

namespace {
// Tag of the "objects in parent" tree associated with `parent`
std::string assn_tag(const art::InputTag &parent, const char *suffix) {
  return art::InputTag(parent.label(), parent.instance() + suffix, parent.process()).encode();
}
} // namespace

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
    dump_simides(p.get<bool>("dump_simides", true)),
    simchannel_tag(p.get<art::InputTag>("simchannel_tag", "tpcrawdecoder:simpleSC"))
{
  need_truth_maps = dump_mctruths || dump_mcparticles || ((dump_tp || dump_ta) && tp_backtracking);

  consumesMany<std::vector<sim::SimChannel>>();

  std::vector<fhicl::ParameterSet> offsets = p.get<std::vector<fhicl::ParameterSet>>("bt_window_offsets");
  std::map<std::string, ViewOffsets> bt_offsets;
  for (const auto &offset : offsets) {
    bt_offsets[offset.get<std::string>("tool_type")] =
        ViewOffsets{offset.get<int>("U"), offset.get<int>("V"), offset.get<int>("X")};
  }
  if (tp_backtracking) {
    tp_bt_ = std::make_unique<TPBacktracker>(std::move(bt_offsets));
  }
}

void dunetrigger::TriggerAnaTree::beginJob() {
  if (dump_mctruths) {

    mctruth_tree = tfs->make<TTree>("mctruths", "mctruths");
    ev_sbuf.make_branches(*mctruth_tree);
    mctruth_buffer.make_branches(*mctruth_tree);


    mcneutrino_tree = tfs->make<TTree>("mcneutrinos", "mcneutrinos");
    ev_sbuf.make_branches(*mcneutrino_tree);
    mcneutrino_buffer.make_branches(*mcneutrino_tree);

  }

  if (dump_mcparticles) {

    mcparticle_tree = tfs->make<TTree>("mcparticles", "mcparticles");
    ev_sbuf.make_branches(*mcparticle_tree);    
    mcparticle_buffer.make_branches(*mcparticle_tree);
  }

  if (dump_simides) {

    simide_tree = tfs->make<TTree>("simides", "simides");
    ev_sbuf.make_branches(*simide_tree);    
    simide_buffer.make_branches(*simide_tree);
  }

  // Summary Trees - always created
  summary_tree = tfs->make<TTree>("event_summary", "event_summary");
  ev_sbuf.make_branches(*summary_tree);
  evsummary_buf.make_branches(*summary_tree);

  simide_summary_tree = tfs->make<TTree>("simide_summary", "simide_summary");
  ev_sbuf.make_branches(*simide_summary_tree);
  simide_summary_buffer.make_branches(*simide_summary_tree);
  simide_tpc_buffer.make_branches(*simide_summary_tree);

  // Save detector settings
  auto const &geo = art::ServiceHandle<geo::Geometry>();

  // Geometry
  info_data["geo"] = {};
  info_data["geo"]["detector"] = geo->DetectorName();

  // Detector Properties
  auto const detProp =
      art::ServiceHandle<detinfo::DetectorPropertiesService const>{}
          ->DataForJob();

  info_data["detector_properties"]["electrons_to_adc"] = detProp.ElectronsToADC();
  info_data["detector_properties"]["electron_lifetime"] = detProp.ElectronLifetime();
  info_data["detector_properties"]["readout_window"] = detProp.ReadOutWindowSize();
  info_data["detector_properties"]["drift_velocity"] = detProp.DriftVelocity();  

  // Backtracking
  if (tp_backtracking) {
    for (const auto &[tool, offsets] : tp_bt_->offsets()) {
      info_data["backtracker"][tool]["offset_U"] = offsets.u;
      info_data["backtracker"][tool]["offset_V"] = offsets.v;
      info_data["backtracker"][tool]["offset_X"] = offsets.w;
    }
  }

  first_event_flag = true;
}

void dunetrigger::TriggerAnaTree::analyze(art::Event const &e) {

  ev_sbuf.reset();
  ev_sbuf->run = e.run();
  ev_sbuf->subrun = e.subRun();
  ev_sbuf->event = e.event();

  evsummary_buf.reset();
  mctruth_buffer.clear();
  mcneutrino_buffer.clear();
  mcparticle_buffer.clear();
  simide_buffer.clear();
  simide_summary_buffer.reset();
  simide_tpc_buffer.clear();
  track_en_sums.clear();
  track_electron_sums.clear();
  simide_tpc_energy_map.clear();
  bt_map.clear();
  tp_tool_type_cache_.clear();
  mctruth_handles.clear();
  trkId_to_truthBlockId.clear();
  truthBlockId_to_generator_name.clear();

  // Clear all TP, TA and TC writers
  for (auto &[tag, ts] : tp_trees) ts.writers.clear();
  for (auto &[tag, ts] : ta_trees) ts.writers.clear();
  for (auto &[tag, ts] : tc_trees) ts.writers.clear();

  // Counters are incremented by the fill* functions below
  evsummary_buf->mctruths_count = 0;
  evsummary_buf->mcparticles_count = 0;
  evsummary_buf->mcneutrinos_count = 0;
  evsummary_buf->simides_count = 0;

  // get a service handle for geometry
  geo::WireReadoutGeom const *geom = &art::ServiceHandle<geo::WireReadout>()->Get();

  if (need_truth_maps) build_truth_maps(e);
  if (dump_mctruths) fillMCTruth(e);
  fillSimChannels(e, geom);
  fillSimIDESummary();
  if (dump_mcparticles) fillMCParticles(e);
  if (dump_tp) fillTPs(e, geom);
  if (dump_ta) fillTAs(e, geom);
  if (dump_tc) fillTCs(e);

  summary_tree->Fill();

  first_event_flag = false;
}

// build_truth_maps
// Reads:    simb::MCTruth collections (all), MCTruth->MCParticle assns from largeant
// Produces: mctruth_handles, trkId_to_truthBlockId, truthBlockId_to_generator_name,
//           info_data["mctruth_blockid_map"], evsummary_buf (mctruths_count, mcneutrinos_count)
// Requires: nothing
void dunetrigger::TriggerAnaTree::build_truth_maps(art::Event const &e) {
  mctruth_handles = e.getMany<std::vector<simb::MCTruth>>();

  int truth_block_counter = 0;

  for (auto const &mctruthHandle : mctruth_handles) {
    // Extract the generator name from the truth handle input label
    std::string generator_name = mctruthHandle.provenance()->inputTag().label();
    // Store generator name for TP backtracking
    truthBlockId_to_generator_name[truth_block_counter] = generator_name;

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
        ++evsummary_buf->mcneutrinos_count;
      }
      evsummary_buf->mctruths_count += truthblock.NParticles();
      truth_block_counter++;
    }
  }

  json j_mctruth_gen_map(truthBlockId_to_generator_name);
  info_data["mctruth_blockid_map"] = j_mctruth_gen_map;
}

// fillMCTruth
// Reads:    mctruth_handles
// Produces: mctruth_buffer/mctruth_tree, mcneutrino_buffer/mcneutrino_tree
// Requires: mctruth_handles (build_truth_maps)
void dunetrigger::TriggerAnaTree::fillMCTruth(art::Event const &e) {
  int truth_block_counter = 0;

  size_t mctruth_collection_size{0};
  for (auto const &mctruthHandle : mctruth_handles) {
    for (size_t i = 0; i < mctruthHandle->size(); i++) {
      const simb::MCTruth &truthblock = *art::Ptr<simb::MCTruth>(mctruthHandle, i);
      mctruth_collection_size += truthblock.NParticles();
    }
  }


  mctruth_buffer.reserve(mctruth_collection_size);

  for (auto const &mctruthHandle : mctruth_handles) {
    // Extract the generator name from the truth handle input label
    std::string generator_name = mctruthHandle.provenance()->inputTag().label();

    for (size_t i = 0; i < mctruthHandle->size(); i++) {
      const simb::MCTruth &truthblock = *art::Ptr<simb::MCTruth>(mctruthHandle, i);
      if (truthblock.NeutrinoSet()) {

        const simb::MCNeutrino &mcneutrino = truthblock.GetNeutrino();
        

        mcneutrino_buffer->block_id = truth_block_counter;
        mcneutrino_buffer->generator_name = generator_name;
        mcneutrino_buffer->nupdg = mcneutrino.Nu().PdgCode();
        mcneutrino_buffer->leptonpdg = mcneutrino.Lepton().PdgCode();
        mcneutrino_buffer->ccnc = mcneutrino.CCNC();
        mcneutrino_buffer->mode = mcneutrino.Mode();
        mcneutrino_buffer->interactionType = mcneutrino.InteractionType();
        mcneutrino_buffer->target = mcneutrino.Target();
        mcneutrino_buffer->hitnuc = mcneutrino.HitNuc();
        mcneutrino_buffer->hitquark = mcneutrino.HitQuark();
        mcneutrino_buffer->w = mcneutrino.W();
        mcneutrino_buffer->x = mcneutrino.X();
        mcneutrino_buffer->y = mcneutrino.Y();
        mcneutrino_buffer->qsqr = mcneutrino.QSqr();
        mcneutrino_buffer->pt = mcneutrino.Pt();
        mcneutrino_buffer->theta = mcneutrino.Theta();
        mcneutrino_buffer.push_back();

      }

      int nparticles = truthblock.NParticles();


      for (int ipart = 0; ipart < nparticles; ipart++) {

        const simb::MCParticle &part = truthblock.GetParticle(ipart);

        mctruth_buffer->block_id = truth_block_counter;
        mctruth_buffer->pdg = part.PdgCode();
        mctruth_buffer->generator_name = generator_name;
        mctruth_buffer->status_code = part.StatusCode();
        mctruth_buffer->process = part.Process();
        mctruth_buffer->truth_track_id = part.TrackId();
        mctruth_buffer->x = part.Vx();
        mctruth_buffer->y = part.Vy();
        mctruth_buffer->z = part.Vz();
        mctruth_buffer->t = part.T();
        mctruth_buffer->px = part.Px();
        mctruth_buffer->py = part.Py();
        mctruth_buffer->pz = part.Pz();
        mctruth_buffer->p = part.P();
        mctruth_buffer->energy = part.E();
        mctruth_buffer->kinetic_energy = part.E() - part.Mass();

        mctruth_buffer.push_back();
      }
      truth_block_counter++;
    }
  }

  mcneutrino_tree->Fill();

  mctruth_tree->Fill();
}

// fillSimChannels
// Reads:    sim::SimChannel collections matching simchannel_tag
// Produces: track_en_sums, track_electron_sums, simide_tpc_energy_map,
//           bt_map, evsummary_buf (per-ROP sums, simides_count),
//           simide_buffer/simide_tree
// Requires: nothing
void dunetrigger::TriggerAnaTree::fillSimChannels(art::Event const &e, geo::WireReadoutGeom const *geom) {

  auto simchannels_many =
      lar::util::getManyByRegexTag<std::vector<sim::SimChannel>>(e, simchannel_tag);
  if (simchannels_many.empty()) {
    throw cet::exception("TriggerAnaTree")
        << "Found no std::vector<sim::SimChannel> collections matching simchannel_tag \""
        << simchannel_tag.encode() << "\"";
  }

  // TODO: alternative implementation that does not rely on `wcls_main.structs.process_apa_index`
  // Get the number of TPCSets from the wiregeometry
  // Loop from 0 to NTPCSets
  // getValidHandle("simpleSC{i_tpcset}")
  // if doesn't exist -> handle
  // else continue as it is

  for ( auto simchannels : simchannels_many ) {

    std::set<int> tpcset_ids;

    for (const sim::SimChannel &sc : *simchannels) {

      
      ChannelInfo chinfo = get_channel_info_for_channel(geom, sc.Channel());
      // Track what TPC elements are in this collection
      tpcset_ids.insert(chinfo.tpcset_id);
      
      
      sim::SimChannel::TDCIDEs_t const &tdcidemap = sc.TDCIDEMap();

      for (const sim::TDCIDE &tdcide : tdcidemap) {
        for (const sim::IDE& ide : tdcide.second) {

          track_en_sums[ide.trackID] += ide.energy;
          track_electron_sums[ide.trackID] += ide.numElectrons;

          // save visible energy only in collection views (for ROI studies)
          if (chinfo.view == geo::kW) {
            simide_tpc_energy_map[chinfo].energy  += ide.energy;
            simide_tpc_energy_map[chinfo].num_electrons += ide.numElectrons;
          }

          // populate per-plane visible energy counters
          if (chinfo.rop_id == 0) {
            evsummary_buf->tot_visible_energy_rop0 += ide.energy;
            evsummary_buf->tot_numelectrons_rop0 += ide.numElectrons;
          }
          else if (chinfo.rop_id == 1) {
            evsummary_buf->tot_visible_energy_rop1 += ide.energy;
            evsummary_buf->tot_numelectrons_rop1 += ide.numElectrons;
          }
          else if (chinfo.rop_id == 2) {
            evsummary_buf->tot_visible_energy_rop2 += ide.energy;
            evsummary_buf->tot_numelectrons_rop2 += ide.numElectrons;
          }
          else if (chinfo.rop_id == 3) {
            evsummary_buf->tot_visible_energy_rop3 += ide.energy;
            evsummary_buf->tot_numelectrons_rop3 += ide.numElectrons;
          }

          if (dump_simides) {
            simide_buffer->channel = sc.Channel();
            simide_buffer->timestamp = tdcide.first;
            simide_buffer->numelectrons = ide.numElectrons;
            simide_buffer->energy = ide.energy;
            simide_buffer->x = ide.x;
            simide_buffer->y = ide.y;
            simide_buffer->z = ide.z;
            simide_buffer->trackID = ide.trackID;
            simide_buffer->origTrackID = ide.origTrackID;
            simide_buffer->readout_plane_id = chinfo.rop_id;
            simide_buffer->readout_view = chinfo.view;
            simide_buffer->detector_element = chinfo.tpcset_id;
            simide_buffer.push_back();
            ++evsummary_buf->simides_count;
          }
        }
      }
    }

    // Each TPCSet must be covered by a single SimChannel collection,
    // otherwise backtracking would silently use only the last one.
    std::vector<int> overlap;
    boost::set_intersection(tpcset_ids, bt_map | boost::adaptors::map_keys,
                            std::back_inserter(overlap));
    if (!overlap.empty()) {
      cet::exception ex("TriggerAnaTree");
      ex << "SimChannel collection " << simchannels.provenance()->inputTag().encode()
         << " covers TPCSets already provided by another collection:";
      for (int tpcset_id : overlap) ex << " " << tpcset_id;
      ex << ". Check that simchannel_tag matches only one set of SimChannel collections.";
      throw ex;
    }

    auto mtb = std::make_shared<MiniBackTracker>(simchannels);
    for (int tpcset_id : tpcset_ids) {
      bt_map[tpcset_id] = mtb;
    }
  }

  if (dump_simides) simide_tree->Fill();
}

// fillSimIDESummary
// Reads:    nothing from the event
// Produces: simide_summary_buffer, simide_tpc_buffer/simide_summary_tree
// Requires: simide_tpc_energy_map (fillSimChannels)
void dunetrigger::TriggerAnaTree::fillSimIDESummary() {
  // Fill simide_summary_tree: one row per {rop, tpcset} for collection-view channels
  double total_visible_energy = 0.;
  double total_numelectrons = 0.;
  for (const auto& [chinfo, edep] : simide_tpc_energy_map) {
    total_visible_energy += edep.energy;
    total_numelectrons   += edep.num_electrons;
  }
  for (const auto& [chinfo, edep] : simide_tpc_energy_map) {
    simide_summary_buffer->total_visible_energy = total_visible_energy;
    simide_summary_buffer->total_numelectrons   = total_numelectrons;
    simide_tpc_buffer->readout_plane_id     = chinfo.rop_id;
    simide_tpc_buffer->detector_element     = chinfo.tpcset_id;
    simide_tpc_buffer->energy_per_tpc       = edep.energy;
    simide_tpc_buffer->numelectrons_per_tpc = edep.num_electrons;
    simide_tpc_buffer.push_back();
  }
  simide_summary_tree->Fill();
}

// fillMCParticles
// Reads:    simb::MCParticle collections (all)
// Produces: evsummary_buf (mcparticles_count), mcparticle_buffer/mcparticle_tree
// Requires: track_en_sums, track_electron_sums (fillSimChannels);
//           trkId_to_truthBlockId (build_truth_maps)
void dunetrigger::TriggerAnaTree::fillMCParticles(art::Event const &e) {

  std::vector<art::Handle<std::vector<simb::MCParticle>>> mcparticleHandles =
      e.getMany<std::vector<simb::MCParticle>>();

  for (auto const &mcparticleHandle : mcparticleHandles) {

    std::string generator_name = mcparticleHandle.provenance()->inputTag().label();

    for (const simb::MCParticle &part : *mcparticleHandle) {

      mcparticle_buffer->pdg = part.PdgCode();
      mcparticle_buffer->generator_name = generator_name;
      mcparticle_buffer->status_code = part.StatusCode();
      mcparticle_buffer->g4_track_id = part.TrackId();
      mcparticle_buffer->mother = part.Mother();
      int truth_block_id = -1;
      auto it = trkId_to_truthBlockId.find(part.TrackId());
      if (it != trkId_to_truthBlockId.end()) truth_block_id = it->second;
      mcparticle_buffer->truth_block_id = truth_block_id;
      mcparticle_buffer->x = part.Vx();
      mcparticle_buffer->y = part.Vy();
      mcparticle_buffer->z = part.Vz();
      mcparticle_buffer->t = part.T();
      mcparticle_buffer->end_x = part.EndX();
      mcparticle_buffer->end_y = part.EndY();
      mcparticle_buffer->end_z = part.EndZ();
      mcparticle_buffer->end_t = part.EndT();
      mcparticle_buffer->px = part.Px();
      mcparticle_buffer->py = part.Py();
      mcparticle_buffer->pz = part.Pz();
      mcparticle_buffer->energy = part.E();
      mcparticle_buffer->kinetic_energy = part.E() - part.Mass();
      mcparticle_buffer->edep = track_en_sums.count(part.TrackId()) ? track_en_sums.at(part.TrackId()) : 0;
      mcparticle_buffer->numelectrons =
          track_electron_sums.count(part.TrackId()) ? track_electron_sums.at(part.TrackId()) : 0;
      mcparticle_buffer->shower_edep = track_en_sums.count(-part.TrackId()) ? track_en_sums.at(-part.TrackId()) : 0;
      mcparticle_buffer->shower_numelectrons =
          track_electron_sums.count(-part.TrackId()) ? track_electron_sums.at(-part.TrackId()) : 0;
      mcparticle_buffer->process = part.Process();
      mcparticle_buffer.push_back();
      ++evsummary_buf->mcparticles_count;
    }
  }
  mcparticle_tree->Fill();
}

// fillTPs
// Reads:    TriggerPrimitive collections matching tp_tag_regex
// Produces: TP trees (tp_trees[<tag>]), info_data["tpg"]
//           (first event only)
// Requires: bt_map (fillSimChannels); trkId_to_truthBlockId,
//           truthBlockId_to_generator_name (build_truth_maps) -- only if tp_backtracking
void dunetrigger::TriggerAnaTree::fillTPs(art::Event const &e, geo::WireReadoutGeom const *geom) {
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

    bool is_tpc_tp_collection = (tp_tool_type.find("TPAlgTPC") == 0);
    if (!is_tpc_tp_collection) {
      throw cet::exception("TriggerAnaTree")
        << "TP collection " << tag << " was produced by tool " << tp_tool_type
        << ", which is not a TPC TP algorithm. Only TPAlgTPC* collections are "
           "supported; adjust tp_tag_regex to exclude it.";
    }

    if ( first_event_flag ) {
      info_data["tpg"][tag]["tool"] = tp_tool_type;

      if (is_tpc_tp_collection) {
        info_data["tpg"][tag]["threshold_tpg_plane0"] = tp_params.get<int>("threshold_tpg_plane0");
        info_data["tpg"][tag]["threshold_tpg_plane1"] = tp_params.get<int>("threshold_tpg_plane1");
        info_data["tpg"][tag]["threshold_tpg_plane2"] = tp_params.get<int>("threshold_tpg_plane2");
      }
    }


    auto &ts = get_or_create_tree(tp_trees, "TriggerPrimitives", "Trigger Primitive Trees",
                                  tag, tp_backtracking, false);

    for (const TriggerPrimitive &tp : *tpHandle) {
      auto chinfo = fill_tp_row(ts.writers, tp, geom);

      // TPC TP backtracking
      if (ts.writers.bt) backtrack_tp(ts.writers, chinfo, tp_tool_type, tag);
    }

    ts.tree->Fill();
  }

}

// fillTAs
// Reads:    TriggerActivityData collections matching ta_tag_regex,
//           TA->TriggerPrimitive assns
// Produces: TA trees (ta_trees[<tag>]), inTAs TP trees (tp_trees[<tag>inTAs])
// Requires: bt_map (fillSimChannels); trkId_to_truthBlockId,
//           truthBlockId_to_generator_name (build_truth_maps) -- only if tp_backtracking
void dunetrigger::TriggerAnaTree::fillTAs(art::Event const &e, geo::WireReadoutGeom const *geom) {
  std::vector<art::Handle<std::vector<TriggerActivityData>>> taHandles =
      e.getMany<std::vector<TriggerActivityData>>();

  std::regex ta_regex(this->ta_tag_regex);

  for (auto const &taHandle : taHandles) {

    art::FindManyP<TriggerPrimitive> assns(taHandle, e, taHandle.provenance()->moduleLabel());
    std::string tag = taHandle.provenance()->inputTag().encode();
    if ( !std::regex_match(tag, ta_regex) ) {
      continue;
    }
    auto &ta_ts = get_or_create_tree(ta_trees, "TriggerActivities", "Trigger Activity Trees",
                                     tag, false, false);
    const std::string tpInTaTag = assn_tag(taHandle.provenance()->inputTag(), "inTAs");
    for (size_t i = 0; i < taHandle->size(); i++) {
      const TriggerActivityData &ta = *art::Ptr<TriggerActivityData>(taHandle, i);
      if (assns.isValid()) {
        size_t ta_idx = i;
        std::vector<art::Ptr<TriggerPrimitive>> matched_tps = assns.at(i);


        auto &tp_ts = get_or_create_tree(tp_trees, "TriggerPrimitives", "Trigger Primitive Trees",
                                         tpInTaTag, tp_backtracking, true);

        for (art::Ptr<TriggerPrimitive> tp : matched_tps) {
          auto chinfo = fill_tp_row(tp_ts.writers, *tp, geom);
          if (tp_ts.writers.bt) backtrack_tp(tp_ts.writers, chinfo, tp_tool_type_for(e, tp.id(), tag), tpInTaTag);
          tp_ts.writers.assn->ta_number = ta_idx;
          tp_ts.writers.assn.push_back();
        }
        tp_ts.tree->Fill();
        tp_ts.writers.clear();
      }
      ta_ts.writers.ta.push_back(ta);
    }
    ta_ts.tree->Fill();
  }
}

// fillTCs
// Reads:    TriggerCandidateData collections matching tc_tag_regex,
//           TC->TriggerActivityData assns
// Produces: TC trees (tc_trees[<tag>]), inTCs TA trees (ta_trees[<tag>inTCs])
// Requires: nothing
void dunetrigger::TriggerAnaTree::fillTCs(art::Event const &e) {
  std::vector<art::Handle<std::vector<TriggerCandidateData>>> tcHandles =
      e.getMany<std::vector<TriggerCandidateData>>();

  std::regex tc_regex(this->tc_tag_regex);

  for (auto const &tcHandle : tcHandles) {
    art::FindManyP<TriggerActivityData> assns(tcHandle, e, tcHandle.provenance()->moduleLabel());
    std::string tag = tcHandle.provenance()->inputTag().encode();
    if ( !std::regex_match(tag, tc_regex) ) {
      continue;
    }
    auto &tc_ts = get_or_create_tree(tc_trees, "TriggerCandidates", "Trigger Candidate Trees",
                                     tag, false, false);
    for (size_t i = 0; i < tcHandle->size(); i++) {
      const TriggerCandidateData &tc = *art::Ptr<TriggerCandidateData>(tcHandle, i);
      if (assns.isValid()) {
        art::InputTag tc_input_tag = tcHandle.provenance()->inputTag();
        auto &ta_ts = get_or_create_tree(ta_trees, "TriggerActivities", "Trigger Activity Trees",
                                         assn_tag(tc_input_tag, "inTCs"), false, true);
        std::vector<art::Ptr<TriggerActivityData>> matched_tas = assns.at(i);
        for (art::Ptr<TriggerActivityData> ta : matched_tas) {
          ta_ts.writers.ta.push_back(*ta);
          ta_ts.writers.assn->tc_number = i;
          ta_ts.writers.assn.push_back();
        }
        ta_ts.tree->Fill();
        ta_ts.writers.ta.clear();
        ta_ts.writers.assn.clear();
      }
      tc_ts.writers.tc.push_back(tc);
    }
    tc_ts.tree->Fill();
  }
}

// fill_tp_row
// Fills the TP staging row from `tp` plus its channel info, commits it, and
// returns the channel info for callers that need it (backtracking).
dunetrigger::ChannelInfo dunetrigger::TriggerAnaTree::fill_tp_row(TPWriters &writers,
                                                                  const TriggerPrimitive &tp,
                                                                  geo::WireReadoutGeom const *geom) {
  auto &tpw = writers.tp;
  tpw->from_tp(tp);
  auto chinfo = get_channel_info_for_channel(geom, tp.channel);
  tpw->readout_plane_id = chinfo.rop_id;
  tpw->readout_view = chinfo.view;
  tpw->TPCSetID = chinfo.tpcset_id;
  tpw.push_back();
  return chinfo;
}

// backtrack_tp
// Fills and commits the backtracking row for the TP last committed to w.tp.
void dunetrigger::TriggerAnaTree::backtrack_tp(TPWriters &w, const ChannelInfo &chinfo,
                                               const std::string &tool_type, const std::string &tag) {
  // SimChannel writers are dense, so every simulated TPCSet must be in
  // bt_map: a missing one means simchannel_tag misses some collections.
  auto bt_it = bt_map.find(chinfo.tpcset_id);
  if (bt_it == bt_map.end()) {
    throw cet::exception("TriggerAnaTree")
        << "No SimChannel collection covers TPCSet " << chinfo.tpcset_id
        << " (TP on channel " << w.tp.row.channel << " from " << tag
        << "). Check that simchannel_tag \"" << simchannel_tag.encode()
        << "\" matches the SimChannels of every TPCSet with TPs.";
  }
  auto& mbt = bt_it->second;
  std::vector<sim::IDE> matched_ides = tp_bt_->match_ides(w.tp.row, tool_type, *mbt);
  tp_bt_->fill_row(w.bt.row, matched_ides, trkId_to_truthBlockId, truthBlockId_to_generator_name, *mbt);
  w.bt.push_back();
}

// tp_tool_type_for
// Returns tpalg.tool_type of the module that produced the TP collection `id`,
// looked up through the product's provenance and cached for the event.
const std::string &dunetrigger::TriggerAnaTree::tp_tool_type_for(art::Event const &e, art::ProductID id,
                                                                 const std::string &ta_tag) {
  auto it = tp_tool_type_cache_.find(id);
  if (it != tp_tool_type_cache_.end()) return it->second;

  auto prov = e.getProductProvenance(id);
  if (!prov || !prov->isValid()) {
    throw cet::exception("TriggerAnaTree")
        << "Cannot find the provenance of TP collection " << id
        << " associated to TA collection " << ta_tag
        << "; its TP algorithm tool type cannot be determined.";
  }
  std::string tool_type =
      prov->parameterSet().get<fhicl::ParameterSet>("tpalg").get<std::string>("tool_type");
  if (tool_type.find("TPAlgTPC") != 0) {
    throw cet::exception("TriggerAnaTree")
        << "TP collection " << prov->inputTag().encode() << " associated to TA collection " << ta_tag
        << " was produced by tool " << tool_type
        << ", which is not a TPC TP algorithm. Only TPAlgTPC* collections are supported.";
  }
  return tp_tool_type_cache_.emplace(id, std::move(tool_type)).first->second;
}


void dunetrigger::TriggerAnaTree::endJob() {

  auto n = tfs->make<TNamed>("info", info_data.dump().c_str());
  n->Write();
}


template <typename Writers>
dunetrigger::TriggerAnaTree::TreeSet<Writers> &dunetrigger::TriggerAnaTree::get_or_create_tree(
    std::map<std::string, TreeSet<Writers>> &trees, const std::string &dir_name, const std::string &dir_title,
    const std::string &tag, bool backtracking, bool association) {
  auto it = trees.find(tag);
  if (it != trees.end()) return it->second;

  art::TFileDirectory dir = tfs->mkdir(dir_name, dir_title);
  mf::LogInfo("TriggerAnaTree") << "Creating new TTree for " << tag;

  // Replace ":" with "_" in TTree names so that they can be used in ROOT's
  // intepreter
  std::string tree_name = tag;
  std::replace(tree_name.begin(), tree_name.end(), ':', '_');

  auto &ts = trees[tag];
  ts.tree = dir.make<TTree>(tree_name.c_str(), tree_name.c_str());
  ev_sbuf.make_branches(*ts.tree);
  ts.writers.enable(backtracking, association);
  ts.writers.make_branches(*ts.tree);   // disabled writers add no branches
  return ts;
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

// ---------------------------------------------------------------------------
// Out-of-line method implementations for structs declared in
// TriggerAnaTree_module.hh
// ---------------------------------------------------------------------------

void dunetrigger::TriggerPrimitiveRow::from_tp(const dunedaq::trgdataformats::TriggerPrimitive &tp) {
  version = 2; // temp, since variables below are converted to v2 version while the TP version in TriggerSim is still 1.
               // Go back to "= tp.version" after changing triggeralgs to v5 (and using TriggerPrimitive2.hpp as header)
  flag = 0;
  detid = tp.detid;
  channel = tp.channel;
  samples_over_threshold = tp.time_over_threshold / dunetrigger::TPAlgTPCTool::ADC_SAMPLING_RATE_IN_DTS;
  time_start = tp.time_start;
  samples_to_peak = (tp.time_peak - tp.time_start) / dunetrigger::TPAlgTPCTool::ADC_SAMPLING_RATE_IN_DTS;
  adc_integral = tp.adc_integral;
  adc_peak = tp.adc_peak;
}

DEFINE_ART_MODULE(dunetrigger::TriggerAnaTree)
