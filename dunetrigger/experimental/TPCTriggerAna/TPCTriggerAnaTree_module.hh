#ifndef DUNETRIGGER_EXPERIMENTAL_TPCTRIGGERANATREE_MODULE_HH
#define DUNETRIGGER_EXPERIMENTAL_TPCTRIGGERANATREE_MODULE_HH

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "larcore/Geometry/WireReadout.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"

#include "detdataformats/trigger/TriggerPrimitive.hpp"
#include "dunetrigger/experimental/TPCTriggerObjs/TriggerActivityData.hpp"
#include "dunetrigger/experimental/TPCTriggerObjs/TriggerActivityROPData.hpp"

// Tree-writing helpers and row types shared with TriggerAna
// (EventMetaData, ChannelInfo, TriggerPrimitiveRow, field-name reflection, buffers)
#include "dunetrigger/TriggerAna/FieldNames.hh"
#include "dunetrigger/TriggerAna/ScalarFieldsBuffer.hh"
#include "dunetrigger/TriggerAna/TriggerAnaRows.hh"
#include "dunetrigger/TriggerAna/VectorFieldsBuffer.hh"

#include <TTree.h>

#include <map>
#include <string>

namespace dunetrigger {

// Field names of the experimental TAs, used as branch names (C++17 path; no-op in C++20).
// Must be in namespace dunetrigger, where the FieldNames trait lives.
REGISTER_FIELD_NAMES(dunedaq::trgdataformats::exp::TriggerActivityROPData,
                     version, time_start, time_end, time_peak, time_activity,
                     channel_start, channel_end, channel_peak,
                     adc_integral, adc_peak, detid, type, algorithm);

REGISTER_FIELD_NAMES(dunedaq::trgdataformats::exp::TriggerActivityData,
                     version, time_start, time_end);

namespace experimental {

//-----------------------------------------------------------------------
// Per-event summary stored next to EventMetaData in event_summary
struct TPCEventSummaryData {
  int ta_rop_collections_count = 0; // TriggerActivityROPData collections matching ta_rop_tag_regex in this event
  int ta_rops_count = 0;            // total TriggerActivityROPData objects over those collections
  int ta_collections_count = 0;     // TriggerActivityData collections matching ta_tag_regex in this event
  int tas_count = 0;                // total TriggerActivityData objects over those collections
};

//-----------------------------------------------------------------------
// Per-ROP-TA association info, one element per ROP TA in the ROP TA tree
struct TriggerActivityROPTPCountRow {
  int n_tps = INVALID_NUM; // TPs associated to the ROP TA; INVALID_NUM if the collection has no TA->TP assns
};

//-----------------------------------------------------------------------
// Per-TP association info in the inTAROPs TP trees
struct TriggerPrimitiveTAROPAssociationRow {
  int ta_rop_number = INVALID_NUM; // index of the ROP TA (within its collection and event) the TP belongs to
};

} // namespace experimental

REGISTER_FIELD_NAMES(experimental::TPCEventSummaryData,
                     ta_rop_collections_count,
                     ta_rops_count,
                     ta_collections_count,
                     tas_count);

REGISTER_FIELD_NAMES(experimental::TriggerActivityROPTPCountRow,
                     n_tps);

REGISTER_FIELD_NAMES(experimental::TriggerPrimitiveTAROPAssociationRow,
                     ta_rop_number);

namespace experimental {

//-----------------------------------------------------------------------
// Minimal experimental counterpart of dunetrigger::TriggerAnaTree: writes only the
// experimental TPC per-ROP TriggerActivities (TriggerActivityROPData), their TP
// association counts and the TPs associated to them, the experimental
// TriggerActivities (TriggerActivityData), and an event summary.
class TPCTriggerAnaTree : public art::EDAnalyzer {
public:
  explicit TPCTriggerAnaTree(fhicl::ParameterSet const &p);

  // Plugins should not be copied or assigned.
  TPCTriggerAnaTree(TPCTriggerAnaTree const &) = delete;
  TPCTriggerAnaTree(TPCTriggerAnaTree &&) = delete;
  TPCTriggerAnaTree &operator=(TPCTriggerAnaTree const &) = delete;
  TPCTriggerAnaTree &operator=(TPCTriggerAnaTree &&) = delete;

  // Required functions.
  void beginJob() override;
  void analyze(art::Event const &e) override;

private:
  using TriggerActivityROP = dunedaq::trgdataformats::exp::TriggerActivityROPData;
  using TriggerActivity = dunedaq::trgdataformats::exp::TriggerActivityData;
  using TriggerPrimitive = dunedaq::trgdataformats::TriggerPrimitive;

  using TriggerPrimitiveWriter = VectorFieldsBuffer<TriggerPrimitiveRow>;
  using TriggerPrimitiveTAROPAssociationWriter = VectorFieldsBuffer<TriggerPrimitiveTAROPAssociationRow>;
  using TriggerActivityROPWriter = VectorFieldsBuffer<TriggerActivityROP>;
  using TriggerActivityROPTPCountWriter = VectorFieldsBuffer<TriggerActivityROPTPCountRow>;
  using TriggerActivityWriter = VectorFieldsBuffer<TriggerActivity>;

  // Writers bound to one output tree, per tree kind
  struct TPWriters {
    TriggerPrimitiveWriter                 tp;
    TriggerPrimitiveTAROPAssociationWriter assn;
    void clear() { tp.clear(); assn.clear(); }
    void make_branches(TTree &t) { tp.make_branches(t); assn.make_branches(t); }
  };
  struct TAROPWriters {
    TriggerActivityROPWriter        ta_rop;
    TriggerActivityROPTPCountWriter assn;
    void clear() { ta_rop.clear(); assn.clear(); }
    void make_branches(TTree &t) { ta_rop.make_branches(t); assn.make_branches(t); }
  };
  struct TAWriters {
    TriggerActivityWriter ta;
    void clear() { ta.clear(); }
    void make_branches(TTree &t) { ta.make_branches(t); }
  };
  template <typename Writers>
  struct TreeSet {
    TTree  *tree = nullptr;
    Writers writers;
  };

  art::ServiceHandle<art::TFileService> tfs;

  std::map<std::string, TreeSet<TPWriters>> tp_trees;        // key: <ROP TA tag>inTAROPs
  std::map<std::string, TreeSet<TAROPWriters>> ta_rop_trees; // key: ROP TA tag
  std::map<std::string, TreeSet<TAWriters>> ta_trees;        // key: TA tag

  art::InputTag ta_rop_tag_regex; // selects TriggerActivityROPData collections
  art::InputTag ta_tag_regex;     // selects TriggerActivityData collections

  template <typename Writers>
  TreeSet<Writers> &get_or_create_tree(std::map<std::string, TreeSet<Writers>> &trees,
                                       const std::string &dir_name, const std::string &dir_title,
                                       const std::string &tag);

  static ChannelInfo get_channel_info_for_channel(geo::WireReadoutGeom const *geom, raw::ChannelID_t channel);

  // Per-product fill functions called by analyze(), in this order
  void fillTAROPs(art::Event const &e, geo::WireReadoutGeom const *geom);
  void fillTAs(art::Event const &e);
  void fillEventSummary();

  // Fills and commits one TP row (TP fields plus channel info)
  void fill_tp_row(TPWriters &writers, const TriggerPrimitive &tp, geo::WireReadoutGeom const *geom);

  // Event meta data buffer
  ScalarFieldsBuffer<EventMetaData> ev_sbuf;

  // Event summary
  TTree *summary_tree = nullptr;
  ScalarFieldsBuffer<TPCEventSummaryData> evsummary_buf;
};

} // namespace experimental
} // namespace dunetrigger

#endif // DUNETRIGGER_EXPERIMENTAL_TPCTRIGGERANATREE_MODULE_HH
