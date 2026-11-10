#ifndef DUNETRIGGER_FDTPC_TPCTRIGGERANATREE_MODULE_HH
#define DUNETRIGGER_FDTPC_TPCTRIGGERANATREE_MODULE_HH

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "larcore/Geometry/WireReadout.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"

#include "detdataformats/trigger/TriggerPrimitive.hpp"
#include "dunetrigger/fdtpc/TPCTriggerObjs/TriggerActivityData.hpp"

// Tree-writing helpers and row types shared with TriggerAna
// (EventMetaData, ChannelInfo, TriggerPrimitiveRow, TriggerPrimitiveAssociationRow,
// field-name reflection, buffers)
#include "dunetrigger/TriggerAna/FieldNames.hh"
#include "dunetrigger/TriggerAna/ScalarFieldsBuffer.hh"
#include "dunetrigger/TriggerAna/TriggerAnaRows.hh"
#include "dunetrigger/TriggerAna/VectorFieldsBuffer.hh"

#include <TTree.h>

#include <map>
#include <string>

namespace dunetrigger {

// Field names of the experimental TA, used as branch names (C++17 path; no-op in C++20).
// Must be in namespace dunetrigger, where the FieldNames trait lives.
REGISTER_FIELD_NAMES(dunedaq::trgdataformats::exp::TriggerActivityData,
                     version, time_start, time_end, time_peak, time_activity,
                     channel_start, channel_end, channel_peak,
                     adc_integral, adc_peak, detid, type, algorithm);

namespace fdtpc {

//-----------------------------------------------------------------------
// Per-event summary stored next to EventMetaData in event_summary
struct TPCEventSummaryData {
  int ta_collections_count = 0; // TA collections matching ta_tag_regex in this event
  int tas_count = 0;            // total TAs over those collections
};

//-----------------------------------------------------------------------
// Per-TA association info, one element per TA in the TA tree
struct TriggerActivityTPCountRow {
  int n_tps = INVALID_NUM; // TPs associated to the TA; INVALID_NUM if the collection has no TA->TP assns
};

} // namespace fdtpc

REGISTER_FIELD_NAMES(fdtpc::TPCEventSummaryData,
                     ta_collections_count,
                     tas_count);

REGISTER_FIELD_NAMES(fdtpc::TriggerActivityTPCountRow,
                     n_tps);

namespace fdtpc {

//-----------------------------------------------------------------------
// Minimal fdtpc counterpart of dunetrigger::TriggerAnaTree: writes only the
// experimental TPC TriggerActivities, their TP association counts, the TPs
// associated to them, and an event summary.
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
  using TriggerActivity = dunedaq::trgdataformats::exp::TriggerActivityData;
  using TriggerPrimitive = dunedaq::trgdataformats::TriggerPrimitive;

  using TriggerPrimitiveWriter = VectorFieldsBuffer<TriggerPrimitiveRow>;
  using TriggerPrimitiveAssociationWriter = VectorFieldsBuffer<TriggerPrimitiveAssociationRow>;
  using TriggerActivityWriter = VectorFieldsBuffer<TriggerActivity>;
  using TriggerActivityTPCountWriter = VectorFieldsBuffer<TriggerActivityTPCountRow>;

  // Writers bound to one output tree, per tree kind
  struct TPWriters {
    TriggerPrimitiveWriter            tp;
    TriggerPrimitiveAssociationWriter assn;
    void clear() { tp.clear(); assn.clear(); }
    void make_branches(TTree &t) { tp.make_branches(t); assn.make_branches(t); }
  };
  struct TAWriters {
    TriggerActivityWriter        ta;
    TriggerActivityTPCountWriter assn;
    void clear() { ta.clear(); assn.clear(); }
    void make_branches(TTree &t) { ta.make_branches(t); assn.make_branches(t); }
  };
  template <typename Writers>
  struct TreeSet {
    TTree  *tree = nullptr;
    Writers writers;
  };

  art::ServiceHandle<art::TFileService> tfs;

  std::map<std::string, TreeSet<TPWriters>> tp_trees; // key: <TA tag>inTAs
  std::map<std::string, TreeSet<TAWriters>> ta_trees; // key: TA tag

  art::InputTag ta_tag_regex;

  template <typename Writers>
  TreeSet<Writers> &get_or_create_tree(std::map<std::string, TreeSet<Writers>> &trees,
                                       const std::string &dir_name, const std::string &dir_title,
                                       const std::string &tag);

  static ChannelInfo get_channel_info_for_channel(geo::WireReadoutGeom const *geom, raw::ChannelID_t channel);

  // Per-product fill functions called by analyze(), in this order
  void fillTAs(art::Event const &e, geo::WireReadoutGeom const *geom);
  void fillEventSummary();

  // Fills and commits one TP row (TP fields plus channel info)
  void fill_tp_row(TPWriters &writers, const TriggerPrimitive &tp, geo::WireReadoutGeom const *geom);

  // Event meta data buffer
  ScalarFieldsBuffer<EventMetaData> ev_sbuf;

  // Event summary
  TTree *summary_tree = nullptr;
  ScalarFieldsBuffer<TPCEventSummaryData> evsummary_buf;
};

} // namespace fdtpc
} // namespace dunetrigger

#endif // DUNETRIGGER_FDTPC_TPCTRIGGERANATREE_MODULE_HH
