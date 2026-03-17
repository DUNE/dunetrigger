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

#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/WireReadout.h"
#include "larcoreobj/SimpleTypesAndConstants/readout_types.h"
#include "lardataobj/Simulation/SimChannel.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include <regex>

#include <TDirectory.h>
#include <TFile.h>
#include <TTree.h>


#include "VectorFieldsBuffer.hh"
#include "ScalarFieldsBuffer.hh"
#include "TriggerAnaTree_module.hh"

#include "dunetrigger/TriggerSim/TPAlgTools/TPAlgTPCTool.hh"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

#include <algorithm>
#include <iostream>
#include <map>

#include <nlohmann/json.hpp>
using json = nlohmann::json;

using dunedaq::trgdataformats::TriggerPrimitive;
using dunedaq::trgdataformats::TriggerActivityData;
using dunedaq::trgdataformats::TriggerCandidateData;





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
    simchannel_tag(p.get<std::string>("simchannel_tag", "tpcrawdecoder:simpleSC"))
// More initializers here.
{
  // FIXME: rename `window_offsets` to `bt_window_offsets`
  std::vector<fhicl::ParameterSet> offsets = p.get<std::vector<fhicl::ParameterSet>>("bt_window_offsets");
  for (const auto &offset : offsets) {
    bt_view_offsets[offset.get<std::string>("tool_type")] = {offset.get<int>("U"), offset.get<int>("V"),
                                                             offset.get<int>("X")};
  }
}

void dunetrigger::TriggerAnaTree::beginJob() {
  if (dump_mctruths) {

    mctruth_tree = tfs->make<TTree>("mctruths", "mctruths");
    ev_sbuf.make_branches(*mctruth_tree);
    mctruth_writer.make_branches(*mctruth_tree);


    mcneutrino_tree = tfs->make<TTree>("mcneutrinos", "mcneutrinos");
    ev_sbuf.make_branches(*mcneutrino_tree);
    mcneutrino_writer.make_branches(*mcneutrino_tree);

  }

  if (dump_mcparticles) {

    mcparticle_tree = tfs->make<TTree>("mcparticles", "mcparticles");
    ev_sbuf.make_branches(*mcparticle_tree);    
    mcparticle_writer.make_branches(*mcparticle_tree);
  }

  if (dump_simides) {

    simide_tree = tfs->make<TTree>("simides", "simides");
    ev_sbuf.make_branches(*simide_tree);    
    simide_writer.make_branches(*simide_tree);
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
  

  ev_sbuf.reset();
  ev_sbuf->run = e.run();
  ev_sbuf->subrun = e.subRun();
  ev_sbuf->event = e.event();

  evsummary_buf.reset();
  mctruth_writer.clear();
  mcneutrino_writer.clear();
  mcparticle_writer.clear();
  simide_writer.clear();
  simide_summary_buffer.reset();
  simide_tpc_buffer.clear();
  track_en_sums.clear();
  track_electron_sums.clear();
  simide_tpc_energy_map.clear();

  // Clear all TP writers
  for( auto& [tag, tpw] : tp_writers) {
    std::apply([](auto&... w) { (w.clear(), ...); }, tpw);
  }

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

    size_t mctruth_collection_size{0};
    for (auto const &mctruthHandle : mctruthHandles) {
      for (size_t i = 0; i < mctruthHandle->size(); i++) {
        const simb::MCTruth &truthblock = *art::Ptr<simb::MCTruth>(mctruthHandle, i);
        mctruth_collection_size += truthblock.NParticles();
      }
    }


    mctruth_writer.reserve(mctruth_collection_size);

    for (auto const &mctruthHandle : mctruthHandles) {
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

          const simb::MCNeutrino &mcneutrino = truthblock.GetNeutrino();
          

          mcneutrino_writer->block_id = truth_block_counter;
          mcneutrino_writer->generator_name = generator_name;
          mcneutrino_writer->nupdg = mcneutrino.Nu().PdgCode();
          mcneutrino_writer->leptonpdg = mcneutrino.Lepton().PdgCode();
          mcneutrino_writer->ccnc = mcneutrino.CCNC();
          mcneutrino_writer->mode = mcneutrino.Mode();
          mcneutrino_writer->interactionType = mcneutrino.InteractionType();
          mcneutrino_writer->target = mcneutrino.Target();
          mcneutrino_writer->hitnuc = mcneutrino.HitNuc();
          mcneutrino_writer->hitquark = mcneutrino.HitQuark();
          mcneutrino_writer->w = mcneutrino.W();
          mcneutrino_writer->x = mcneutrino.X();
          mcneutrino_writer->y = mcneutrino.Y();
          mcneutrino_writer->qsqr = mcneutrino.QSqr();
          mcneutrino_writer->pt = mcneutrino.Pt();
          mcneutrino_writer->theta = mcneutrino.Theta();
          mcneutrino_writer.push_back();


          ++mcneutrinos_count;
        }

        int nparticles = truthblock.NParticles();


        for (int ipart = 0; ipart < nparticles; ipart++) {

          const simb::MCParticle &part = truthblock.GetParticle(ipart);

          mctruth_writer->block_id = truth_block_counter;
          mctruth_writer->pdg = part.PdgCode();
          mctruth_writer->generator_name = generator_name;
          mctruth_writer->status_code = part.StatusCode();
          mctruth_writer->process = part.Process();
          mctruth_writer->truth_track_id = part.TrackId();
          mctruth_writer->x = part.Vx();
          mctruth_writer->y = part.Vy();
          mctruth_writer->z = part.Vz();
          mctruth_writer->t = part.T();
          mctruth_writer->px = part.Px();
          mctruth_writer->py = part.Py();
          mctruth_writer->pz = part.Pz();
          mctruth_writer->p = part.P();
          mctruth_writer->energy = part.E();
          mctruth_writer->kinetic_energy = part.E() - part.Mass();

          mctruth_writer.push_back();

          ++mctruths_count;
        }
        truth_block_counter++;
      }
    }

    mcneutrino_tree->Fill();

    mctruth_tree->Fill();

    json j_mctruth_gen_map(truthBlockId_to_generator_name);
    info_data["mctruth_blockid_map"] = j_mctruth_gen_map;

  }

  {
    auto simchannels = e.getValidHandle<std::vector<sim::SimChannel>>(simchannel_tag);

    for (const sim::SimChannel &sc : *simchannels) {

      ChannelInfo chinfo = get_channel_info_for_channel(geom, sc.Channel());
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
            simide_writer->channel = sc.Channel();
            simide_writer->timestamp = tdcide.first;
            simide_writer->numelectrons = ide.numElectrons;
            simide_writer->energy = ide.energy;
            simide_writer->x = ide.x;
            simide_writer->y = ide.y;
            simide_writer->z = ide.z;
            simide_writer->trackID = ide.trackID;
            simide_writer->origTrackID = ide.origTrackID;
            simide_writer->readout_plane_id = chinfo.rop_id;
            simide_writer->readout_view = chinfo.view;
            simide_writer->detector_element = chinfo.tpcset_id;
            simide_writer.push_back();
            ++simides_count;
          }
        }
      }
    }

    if (dump_simides) simide_tree->Fill();
  }

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

  if (dump_mcparticles) {

    std::vector<art::Handle<std::vector<simb::MCParticle>>> mcparticleHandles =
        e.getMany<std::vector<simb::MCParticle>>();

    for (auto const &mcparticleHandle : mcparticleHandles) {

      std::string generator_name = mcparticleHandle.provenance()->inputTag().label();

      for (const simb::MCParticle &part : *mcparticleHandle) {

        mcparticle_writer->pdg = part.PdgCode();
        mcparticle_writer->generator_name = generator_name;
        mcparticle_writer->status_code = part.StatusCode();
        mcparticle_writer->g4_track_id = part.TrackId();
        mcparticle_writer->mother = part.Mother();
        mcparticle_writer->truth_block_id = dump_mctruths ? trkId_to_truthBlockId.at(part.TrackId()) : -1;
        mcparticle_writer->x = part.Vx();
        mcparticle_writer->y = part.Vy();
        mcparticle_writer->z = part.Vz();
        mcparticle_writer->t = part.T();
        mcparticle_writer->end_x = part.EndX();
        mcparticle_writer->end_y = part.EndY();
        mcparticle_writer->end_z = part.EndZ();
        mcparticle_writer->end_t = part.EndT();
        mcparticle_writer->px = part.Px();
        mcparticle_writer->py = part.Py();
        mcparticle_writer->pz = part.Pz();
        mcparticle_writer->energy = part.E();
        mcparticle_writer->kinetic_energy = part.E() - part.Mass();
        mcparticle_writer->edep = track_en_sums.count(part.TrackId()) ? track_en_sums.at(part.TrackId()) : 0;
        mcparticle_writer->numelectrons =
            track_electron_sums.count(part.TrackId()) ? track_electron_sums.at(part.TrackId()) : 0;
        mcparticle_writer->shower_edep = track_en_sums.count(-part.TrackId()) ? track_en_sums.at(-part.TrackId()) : 0;
        mcparticle_writer->shower_numelectrons =
            track_electron_sums.count(-part.TrackId()) ? track_electron_sums.at(-part.TrackId()) : 0;
        mcparticle_writer->process = part.Process();
        mcparticle_writer.push_back();
        ++mcparticles_count;
      }
    }
    mcparticle_tree->Fill();
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

      bool is_tpc_tp_collection = (tp_tool_type.find("TPAlgTPC") == 0);
      // bool is_pds_tp_collection = (tp_tool_type.find("TPAlgPDS") == 0);

      if ( first_event_flag ) {
        info_data["tpg"][tag]["tool"] = tp_tool_type;

        if (is_tpc_tp_collection) {
          info_data["tpg"][tag]["threshold_tpg_plane0"] = tp_params.get<int>("threshold_tpg_plane0");
          info_data["tpg"][tag]["threshold_tpg_plane1"] = tp_params.get<int>("threshold_tpg_plane1");
          info_data["tpg"][tag]["threshold_tpg_plane2"] = tp_params.get<int>("threshold_tpg_plane2");
        }
      }


      std::string map_tag = "tp/" + tag;

      make_tp_tree_if_needed(tag);

      TTree *tp_tree = tree_map[map_tag];

      auto& [tp_writer, tpbt_writer, tpass_writer] = tp_writers[map_tag];

      for (const TriggerPrimitive &tp : *tpHandle) {
        tp_writer->from_tp(tp);
        auto chinfo = get_channel_info_for_channel(geom, tp.channel);
        tp_writer->readout_plane_id = chinfo.rop_id;
        tp_writer->readout_view = chinfo.view;
        tp_writer->TPCSetID = chinfo.tpcset_id;
        tp_writer.push_back();

        // TPC TP backtracking
        if (tpbt_writer and is_tpc_tp_collection) {
          std::vector<sim::IDE> matched_ides = match_simides_to_tps(tp_writer.row, tp_tool_type);
          tpbt_writer->populate_backtracking_info(matched_ides, trkId_to_truthBlockId, truthBlockId_to_generator_name);
          tpbt_writer.push_back();
        }
      }

      tp_tree->Fill();
    }

  }

  evsummary_buf->mctruths_count = mctruths_count;
  evsummary_buf->mcparticles_count = mcparticles_count;
  evsummary_buf->mcneutrinos_count = mcneutrinos_count;
  evsummary_buf->simides_count = simides_count;

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
          fAssnIdx = i;
          std::vector<art::Ptr<TriggerPrimitive>> matched_tps = assns.at(i);


          std::string map_tpInTaTag = "tp/" + tpInTaTag;
          make_tp_tree_if_needed(tpInTaTag, true);
          TTree *tp_tree = tree_map[map_tpInTaTag];
          auto& [tp_writer, tpbt_writer, tpass_writer] = tp_writers[map_tpInTaTag];

          for (art::Ptr<TriggerPrimitive> tp : matched_tps) {
            tp_writer->from_tp(*tp);
            auto chinfo = get_channel_info_for_channel(geom, tp->channel);
            tp_writer->readout_plane_id = chinfo.rop_id;
            tp_writer->readout_view = chinfo.view;
            tp_writer->TPCSetID = chinfo.tpcset_id;
            tp_writer.push_back();
            if (tpbt_writer) tpbt_writer.push_back(); // push default (INVALID) row -- backtracking not computed for association TPs
            tpass_writer->ta_number = fAssnIdx;
            tpass_writer.push_back();
          }
          tp_tree->Fill();
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

  summary_tree->Fill();

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
    TTree* tree = tp_dir.make<TTree>(tree_name.c_str(), tree_name.c_str());
    tree_map[map_tag] = tree;


    ev_sbuf.make_branches(*tree);

    auto& [tpw, tpbtw, tpassw] = tp_writers[map_tag];
    tpbtw.enable(tp_backtracking);
    tpassw.enable(assn);
    tpw.make_branches(*tree);
    tpbtw.make_branches(*tree);   // no-op if disabled
    tpassw.make_branches(*tree);

    // }
    // auto& [curr_tp_writer, curr_tpbt_writer] = it->second;
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

std::vector<sim::IDE> dunetrigger::TriggerAnaTree::match_simides_to_tps(const TriggerPrimitiveRow &tp,
                                                                        const std::string &tool_type) const {

  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  auto it = bt_view_offsets.find(tool_type);
  if (it == bt_view_offsets.end()) {
    std::cout << "No offsets found for tool type " << tool_type << ", using 0,0,0" << std::endl;
  }
  const std::array<int, 3> &offsets = it != bt_view_offsets.end() ? it->second : std::array<int, 3>{0, 0, 0};
  int offset = 0;
  switch (tp.readout_view) {
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

void dunetrigger::TriggerPrimitiveBacktrackingRow::populate_backtracking_info(
    const std::vector<sim::IDE> &ides,
    const std::unordered_map<int, int> &trkid_to_truth_block,
    const std::unordered_map<int, std::string> &truth_id_to_gen) {
  bt_primary_track_id = INVALID;
  bt_primary_track_numelectron_frac = INVALID;
  bt_primary_track_energy_frac = INVALID;
  bt_edep = 0.;
  bt_numelectrons = 0.;
  bt_x = INVALID;
  bt_y = INVALID;
  bt_z = INVALID;
  bt_primary_x = INVALID;
  bt_primary_y = INVALID;
  bt_primary_z = INVALID;
  bt_truth_block_id = INVALID;
  bt_generator_name.clear();

  if (ides.empty()) {
    return;
  }

  std::vector<sim::IDE> bt_ides;
  bt_ides.reserve(ides.capacity());
  std::copy_if(ides.begin(),
               ides.end(),
               std::back_inserter(bt_ides),
               [](const sim::IDE &ide) { return ide.trackID != 0; });

  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;

  std::map<int, double> track_numelectrons;
  std::map<int, double> track_energies;

  if (bt_ides.empty()) {
    std::cout << "Empty IDEs set!" << std::endl;
    return;
  }

  for (const sim::IDE &ide : bt_ides) {
    int mc_track_id = pi_serv->TrackIdToParticle_P(ide.trackID)->TrackId();
    track_numelectrons[mc_track_id] += ide.numElectrons;
    track_energies[mc_track_id] += ide.energy;
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

  bt_truth_block_id = trkid_to_truth_block.at(bt_primary_track_id);
  bt_generator_name = truth_id_to_gen.at(bt_truth_block_id);
}

// ---------------------------------------------------------------------------
// C++17 field name registration.
// To be removed when LArSoft switches to C++20 standard.
// ---------------------------------------------------------------------------

REGISTER_FIELD_NAMES(dunetrigger::MCTruthRow,
                         pdg,
                         process,
                         status_code,
                         block_id,
                         truth_track_id,
                         generator_name,
                         x,
                         y,
                         z,
                         t,
                         px,
                         py,
                         pz,
                         p,
                         energy,
                         kinetic_energy)

REGISTER_FIELD_NAMES(dunetrigger::MCNeutrinoRow,
                         block_id,
                         generator_name,
                         nupdg,
                         leptonpdg,
                         ccnc,
                         mode,
                         interactionType,
                         target,
                         hitnuc,
                         hitquark,
                         w,
                         x,
                         y,
                         qsqr,
                         pt,
                         theta)

REGISTER_FIELD_NAMES(dunetrigger::MCParticleRow,
                         pdg,
                         generator_name,
                         status_code,
                         g4_track_id,
                         mother,
                         truth_block_id,
                         x,
                         y,
                         z,
                         t,
                         end_x,
                         end_y,
                         end_z,
                         end_t,
                         px,
                         py,
                         pz,
                         energy,
                         kinetic_energy,
                         edep,
                         numelectrons,
                         shower_edep,
                         shower_numelectrons,
                         process)

REGISTER_FIELD_NAMES(dunetrigger::SimIDERow,
                         channel,
                         timestamp,
                         numelectrons,
                         energy,
                         x,
                         y,
                         z,
                         trackID,
                         origTrackID,
                         readout_plane_id,
                         readout_view,
                         detector_element)

REGISTER_FIELD_NAMES(dunetrigger::SimIDESummaryRow,
                         total_visible_energy,
                         total_numelectrons)

REGISTER_FIELD_NAMES(dunetrigger::SimIDETPCRow,
                         readout_plane_id,
                         detector_element,
                         energy_per_tpc,
                         numelectrons_per_tpc)

REGISTER_FIELD_NAMES(dunetrigger::TriggerPrimitiveRow,
                         version,
                         flag,
                         detid,
                         channel,
                         samples_over_threshold,
                         time_start,
                         samples_to_peak,
                         adc_integral,
                         adc_peak,
                         readout_plane_id,
                         readout_view,
                         TPCSetID)

REGISTER_FIELD_NAMES(dunetrigger::TriggerPrimitiveBacktrackingRow,
                         bt_primary_track_id,
                         bt_primary_track_numelectron_frac,
                         bt_primary_track_energy_frac,
                         bt_edep,
                         bt_numelectrons,
                         bt_x,
                         bt_y,
                         bt_z,
                         bt_primary_x,
                         bt_primary_y,
                         bt_primary_z,
                         bt_truth_block_id,
                         bt_generator_name)

REGISTER_FIELD_NAMES(dunetrigger::TriggerPrimitiveAssociationRow,
                         ta_number)

REGISTER_FIELD_NAMES(dunetrigger::EventMetaData,
                            event,
                            run,
                            subrun)

REGISTER_FIELD_NAMES(dunetrigger::EventSummaryData,
                            mctruths_count,
                            mcparticles_count,
                            mcneutrinos_count,
                            simides_count,
                            tot_visible_energy_rop0,
                            tot_visible_energy_rop1,
                            tot_visible_energy_rop2,
                            tot_visible_energy_rop3,
                            tot_numelectrons_rop0,
                            tot_numelectrons_rop1,
                            tot_numelectrons_rop2,
                            tot_numelectrons_rop3)

DEFINE_ART_MODULE(dunetrigger::TriggerAnaTree)
