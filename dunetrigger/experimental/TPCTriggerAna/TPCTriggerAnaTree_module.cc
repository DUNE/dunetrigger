////////////////////////////////////////////////////////////////////////
// Class:       TPCTriggerAnaTree
// Plugin Type: analyzer
// File:        TPCTriggerAnaTree_module.cc
//
// Minimal experimental counterpart of TriggerAna/TriggerAnaTree: writes only the
// experimental TPC per-ROP TriggerActivities (dunedaq::trgdataformats::exp::TriggerActivityROPData),
// the TPs associated to them, the experimental TriggerActivities
// (dunedaq::trgdataformats::exp::TriggerActivityData), and an event summary tree.
// No truth, SimChannel, standalone TP or TC information.
//
// Output (TFileService), ':' in tags replaced by '_' in tree names:
//   event_summary                       one entry per event: run/subrun/event + ROP TA and TA counts
//   TriggerActivityROPs/<tag>           one tree per TriggerActivityROPData collection matching
//                                       ta_rop_tag_regex; one entry per event, one vector element
//                                       per ROP TA (ROP TA fields + n_tps)
//   TriggerPrimitives/<label>_<instance>inTAROPs_<process>
//                                       TPs associated to the ROP TAs of that collection;
//                                       one entry per ROP TA, one vector element per TP
//                                       (TP fields + channel info + ta_rop_number = ROP TA index)
//   TriggerActivities/<tag>             one tree per TriggerActivityData collection matching
//                                       ta_tag_regex; one entry per event, one vector
//                                       element per TA (version, time_start, time_end)
////////////////////////////////////////////////////////////////////////

#include "TPCTriggerAnaTree_module.hh"

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileDirectory.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "lardata/ArtDataHelper/GetManyByRegexTag.h"
#include "larcoreobj/SimpleTypesAndConstants/readout_types.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <algorithm>
#include <vector>

namespace {
// Tag of the "objects in parent" tree associated with `parent`
std::string assn_tag(const art::InputTag &parent, const char *suffix) {
  return art::InputTag(parent.label(), parent.instance() + suffix, parent.process()).encode();
}
} // namespace

dunetrigger::experimental::TPCTriggerAnaTree::TPCTriggerAnaTree(fhicl::ParameterSet const &p)
    : EDAnalyzer{p},
    ta_rop_tag_regex(p.get<std::string>("ta_rop_tag_regex", ".*")),
    ta_tag_regex(p.get<std::string>("ta_tag_regex", ".*"))
{
  // Products read by analyze()
  consumesMany<std::vector<TriggerActivityROP>>();
  mayConsumeMany<art::Assns<TriggerActivityROP, TriggerPrimitive>>();   // FindManyP ROP TA->TP
  consumesMany<std::vector<TriggerActivity>>();
}

void dunetrigger::experimental::TPCTriggerAnaTree::beginJob() {
  // Summary tree - always created
  summary_tree = tfs->make<TTree>("event_summary", "event_summary");
  ev_sbuf.make_branches(*summary_tree);
  evsummary_buf.make_branches(*summary_tree);
}

void dunetrigger::experimental::TPCTriggerAnaTree::analyze(art::Event const &e) {

  ev_sbuf.reset();
  ev_sbuf->run = e.run();
  ev_sbuf->subrun = e.subRun();
  ev_sbuf->event = e.event();

  evsummary_buf.reset();

  // Clear all TP, ROP TA and TA writers
  for (auto &[tag, ts] : tp_trees) ts.writers.clear();
  for (auto &[tag, ts] : ta_rop_trees) ts.writers.clear();
  for (auto &[tag, ts] : ta_trees) ts.writers.clear();

  // get a service handle for geometry
  geo::WireReadoutGeom const *geom = &art::ServiceHandle<geo::WireReadout>()->Get();

  fillTAROPs(e, geom);
  fillTAs(e);
  fillEventSummary();
}

// fillTAROPs
// Reads:    exp::TriggerActivityROPData collections matching ta_rop_tag_regex,
//           ROP TA->TriggerPrimitive assns from the same module label
// Produces: ROP TA trees (ta_rop_trees[<tag>]), inTAROPs TP trees (tp_trees[<tag>inTAROPs]),
//           evsummary_buf (ta_rop_collections_count, ta_rops_count)
// Requires: nothing
// Note: trees whose collection is absent from an event are not filled for that event.
//       inTAROPs trees are created and filled only for collections with ROP TA->TP assns.
void dunetrigger::experimental::TPCTriggerAnaTree::fillTAROPs(art::Event const &e, geo::WireReadoutGeom const *geom) {
  auto taRopHandles = lar::util::getManyByRegexTag<std::vector<TriggerActivityROP>>(e, ta_rop_tag_regex);

  for (auto const &taRopHandle : taRopHandles) {

    art::FindManyP<TriggerPrimitive> assns(taRopHandle, e, taRopHandle.provenance()->moduleLabel());
    const std::string tag = taRopHandle.provenance()->inputTag().encode();
    auto &ta_rop_ts = get_or_create_tree(ta_rop_trees, "TriggerActivityROPs", "Trigger Activity ROP Trees", tag);
    const std::string tpInTaRopTag = assn_tag(taRopHandle.provenance()->inputTag(), "inTAROPs");

    ta_rop_ts.writers.ta_rop.reserve(taRopHandle->size());
    ta_rop_ts.writers.assn.reserve(taRopHandle->size());
    for (size_t i = 0; i < taRopHandle->size(); i++) {
      if (assns.isValid()) {
        std::vector<art::Ptr<TriggerPrimitive>> matched_tps = assns.at(i);

        auto &tp_ts = get_or_create_tree(tp_trees, "TriggerPrimitives", "Trigger Primitive Trees", tpInTaRopTag);
        for (art::Ptr<TriggerPrimitive> tp : matched_tps) {
          fill_tp_row(tp_ts.writers, *tp, geom);
          tp_ts.writers.assn->ta_rop_number = i;
          tp_ts.writers.assn.push_back();
        }
        tp_ts.tree->Fill();
        tp_ts.writers.clear();

        ta_rop_ts.writers.assn->n_tps = static_cast<int>(matched_tps.size());
      } else {
        ta_rop_ts.writers.assn->n_tps = INVALID_NUM;
      }
      ta_rop_ts.writers.assn.push_back();
      ta_rop_ts.writers.ta_rop.push_back((*taRopHandle)[i]);
    }
    ta_rop_ts.tree->Fill();

    ++evsummary_buf->ta_rop_collections_count;
    evsummary_buf->ta_rops_count += static_cast<int>(taRopHandle->size());
  }
}

// fillTAs
// Reads:    exp::TriggerActivityData collections matching ta_tag_regex
// Produces: TA trees (ta_trees[<tag>]),
//           evsummary_buf (ta_collections_count, tas_count)
// Requires: nothing
// Note: trees whose collection is absent from an event are not filled for that event.
void dunetrigger::experimental::TPCTriggerAnaTree::fillTAs(art::Event const &e) {
  auto taHandles = lar::util::getManyByRegexTag<std::vector<TriggerActivity>>(e, ta_tag_regex);

  for (auto const &taHandle : taHandles) {
    const std::string tag = taHandle.provenance()->inputTag().encode();
    auto &ta_ts = get_or_create_tree(ta_trees, "TriggerActivities", "Trigger Activity Trees", tag);

    ta_ts.writers.ta.reserve(taHandle->size());
    for (auto const &ta : *taHandle) {
      ta_ts.writers.ta.push_back(ta);
    }
    ta_ts.tree->Fill();

    ++evsummary_buf->ta_collections_count;
    evsummary_buf->tas_count += static_cast<int>(taHandle->size());
  }
}

// fillEventSummary
// Reads:    nothing from the event
// Produces: summary_tree entry
// Requires: evsummary_buf counters (fillTAROPs, fillTAs)
void dunetrigger::experimental::TPCTriggerAnaTree::fillEventSummary() {
  summary_tree->Fill();
}

// fill_tp_row
// Fills the TP staging row from `tp` plus its channel info and commits it.
// The TP conversion is the shared TriggerPrimitiveRow::from_tp (dunetrigger::TriggerAna library).
void dunetrigger::experimental::TPCTriggerAnaTree::fill_tp_row(TPWriters &writers,
                                                        const TriggerPrimitive &tp,
                                                        geo::WireReadoutGeom const *geom) {
  auto &tpw = writers.tp;
  tpw->from_tp(tp);
  auto chinfo = get_channel_info_for_channel(geom, tp.channel);
  tpw->readout_plane_id = chinfo.rop_id;
  tpw->readout_view = chinfo.view;
  tpw->tpcset_id = chinfo.tpcset_id;
  tpw.push_back();
}

template <typename Writers>
dunetrigger::experimental::TPCTriggerAnaTree::TreeSet<Writers> &dunetrigger::experimental::TPCTriggerAnaTree::get_or_create_tree(
    std::map<std::string, TreeSet<Writers>> &trees, const std::string &dir_name, const std::string &dir_title,
    const std::string &tag) {
  auto it = trees.find(tag);
  if (it != trees.end()) return it->second;

  art::TFileDirectory dir = tfs->mkdir(dir_name, dir_title);
  mf::LogInfo("TPCTriggerAnaTree") << "Creating new TTree for " << tag;

  // Replace ":" with "_" in TTree names so that they can be used in ROOT's
  // intepreter
  std::string tree_name = tag;
  std::replace(tree_name.begin(), tree_name.end(), ':', '_');

  auto &ts = trees[tag];
  ts.tree = dir.make<TTree>(tree_name.c_str(), tree_name.c_str());
  ev_sbuf.make_branches(*ts.tree);
  ts.writers.make_branches(*ts.tree);
  return ts;
}

dunetrigger::ChannelInfo dunetrigger::experimental::TPCTriggerAnaTree::get_channel_info_for_channel(
    geo::WireReadoutGeom const *geom, raw::ChannelID_t channel) {
  readout::ROPID rop = geom->ChannelToROP(channel);
  ChannelInfo result;
  result.rop_id = rop.ROP;
  result.tpcset_id = rop.asTPCsetID().TPCset;
  result.view = geom->View(rop);
  return result;
}

DEFINE_ART_MODULE(dunetrigger::experimental::TPCTriggerAnaTree)
