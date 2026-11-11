////////////////////////////////////////////////////////////////////////
// Class:       TPCTriggerAnaTree
// Plugin Type: analyzer
// File:        TPCTriggerAnaTree_module.cc
//
// Minimal experimental counterpart of TriggerAna/TriggerAnaTree: writes only the
// experimental TPC TriggerActivities (dunedaq::trgdataformats::exp::TriggerActivityData),
// the TPs associated to them, and an event summary tree.
// No truth, SimChannel, standalone TP or TC information.
//
// Output (TFileService), ':' in tags replaced by '_' in tree names:
//   event_summary                       one entry per event: run/subrun/event + TA counts
//   TriggerActivities/<tag>             one tree per TA collection matching ta_tag_regex;
//                                       one entry per event, one vector element per TA
//                                       (TA fields + n_tps)
//   TriggerPrimitives/<label>_<instance>inTAs_<process>
//                                       TPs associated to the TAs of that collection;
//                                       one entry per TA, one vector element per TP
//                                       (TP fields + channel info + ta_number = TA index)
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
    ta_tag_regex(p.get<std::string>("ta_tag_regex", ".*"))
{
  // Products read by analyze()
  consumesMany<std::vector<TriggerActivity>>();
  mayConsumeMany<art::Assns<TriggerActivity, TriggerPrimitive>>();   // FindManyP TA->TP
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

  // Clear all TP and TA writers
  for (auto &[tag, ts] : tp_trees) ts.writers.clear();
  for (auto &[tag, ts] : ta_trees) ts.writers.clear();

  // get a service handle for geometry
  geo::WireReadoutGeom const *geom = &art::ServiceHandle<geo::WireReadout>()->Get();

  fillTAs(e, geom);
  fillEventSummary();
}

// fillTAs
// Reads:    exp::TriggerActivityData collections matching ta_tag_regex,
//           TA->TriggerPrimitive assns from the same module label
// Produces: TA trees (ta_trees[<tag>]), inTAs TP trees (tp_trees[<tag>inTAs]),
//           evsummary_buf (ta_collections_count, tas_count)
// Requires: nothing
// Note: trees whose collection is absent from an event are not filled for that event.
//       inTAs trees are created and filled only for collections with TA->TP assns.
void dunetrigger::experimental::TPCTriggerAnaTree::fillTAs(art::Event const &e, geo::WireReadoutGeom const *geom) {
  auto taHandles = lar::util::getManyByRegexTag<std::vector<TriggerActivity>>(e, ta_tag_regex);

  for (auto const &taHandle : taHandles) {

    art::FindManyP<TriggerPrimitive> assns(taHandle, e, taHandle.provenance()->moduleLabel());
    const std::string tag = taHandle.provenance()->inputTag().encode();
    auto &ta_ts = get_or_create_tree(ta_trees, "TriggerActivities", "Trigger Activity Trees", tag);
    const std::string tpInTaTag = assn_tag(taHandle.provenance()->inputTag(), "inTAs");

    ta_ts.writers.ta.reserve(taHandle->size());
    ta_ts.writers.assn.reserve(taHandle->size());
    for (size_t i = 0; i < taHandle->size(); i++) {
      if (assns.isValid()) {
        std::vector<art::Ptr<TriggerPrimitive>> matched_tps = assns.at(i);

        auto &tp_ts = get_or_create_tree(tp_trees, "TriggerPrimitives", "Trigger Primitive Trees", tpInTaTag);
        for (art::Ptr<TriggerPrimitive> tp : matched_tps) {
          fill_tp_row(tp_ts.writers, *tp, geom);
          tp_ts.writers.assn->ta_number = i;
          tp_ts.writers.assn.push_back();
        }
        tp_ts.tree->Fill();
        tp_ts.writers.clear();

        ta_ts.writers.assn->n_tps = static_cast<int>(matched_tps.size());
      } else {
        ta_ts.writers.assn->n_tps = INVALID_NUM;
      }
      ta_ts.writers.assn.push_back();
      ta_ts.writers.ta.push_back((*taHandle)[i]);
    }
    ta_ts.tree->Fill();

    ++evsummary_buf->ta_collections_count;
    evsummary_buf->tas_count += static_cast<int>(taHandle->size());
  }
}

// fillEventSummary
// Reads:    nothing from the event
// Produces: summary_tree entry
// Requires: evsummary_buf counters (fillTAs)
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
