#ifndef DUNETRIGGER_TRIGGERANATREE_MODULE_HH
#define DUNETRIGGER_TRIGGERANATREE_MODULE_HH

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Provenance/ProductID.h"
#include "detdataformats/trigger/TriggerActivityData.hpp"
#include "detdataformats/trigger/TriggerCandidateData.hpp"
#include "detdataformats/trigger/TriggerPrimitive.hpp"
#include "fhiclcpp/ParameterSet.h"
#include "larcore/Geometry/WireReadout.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

#include "FieldNames.hh"
#include "TriggerAnaRows.hh"
#include "ScalarFieldsBuffer.hh"
#include "VectorFieldsBuffer.hh"

#include <TTree.h>

#include <nlohmann/json.hpp>

#include <array>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

namespace dunetrigger {

// Sample-window offsets applied before matching a TP to the SimChannel IDEs,
// one per readout view, indexed by view number (geo::View_t value). DUNE TPCs
// have three readout views: 0 = U, 1 = V, 2 = W/Z (collection); all three are
// required in the configuration.
//
// NOTE (future improvement, not a bug): the right offsets depend on the TP
// algorithm *and its configuration* (thresholds, running-sum window, ...), i.e.
// on the TP collection, not only on the tool type. At the moment one set of
// per-view offsets is configured per tool_type (fcl `bt_window_offsets`) and applied to
// every TP collection produced with that tool, whatever its parameters. Keying
// the offsets by TP collection, or deriving them from the algorithm parameters
// recorded in the producer's provenance, would remove this approximation.
struct ViewOffsets {
  static constexpr int kNumViews = 3;
  std::array<int, kNumViews> by_view{};        // view number -> sample offset
  int for_view(int view) const;                // 0 if view is outside [0, kNumViews)
};

// Matches TPs to SimChannel IDEs and to MC truth. Owns a per-event, merged
// SimChannel index (all matched collections, sorted by channel). The index is
// a stopgap for larsim < v10_27_00, whose cheat::BackTracker cannot load
// several SimChannel collections.
class TPBacktracker {
public:
  explicit TPBacktracker(std::map<std::string, ViewOffsets> offsets);

  // Builds the per-event SimChannel index from all the given collections.
  void prep_event(std::vector<art::Handle<std::vector<sim::SimChannel>>> const& simchannels);
  void clear_event();

  // Returns IDEs on the TP's channel inside the TP's time window shifted by the
  // tool-specific view offset. Empty if no SimChannel exists for the channel.
  std::vector<sim::IDE> match_ides(const TriggerPrimitiveRow& tp,
                                   const std::string& tool_type) const;

  // Fills `row` from `ides`. Leaves INVALID_NUM defaults if `ides` is empty or
  // no IDE could be resolved to an MCParticle.
  void fill_row(TriggerPrimitiveBacktrackingRow& row,
                const std::vector<sim::IDE>& ides,
                const std::unordered_map<int,int>& trkid_to_truth_block,
                const std::unordered_map<int,std::string>& truth_id_to_gen) const;

  const std::map<std::string, ViewOffsets>& offsets() const { return offsets_; }

private:
  // The only two places that touch SimChannels. Step 9b replaces their bodies
  // with cheat::BackTrackerService calls (provider()->FindSimChannelPtr / SimIDEsToXYZ).
  art::Ptr<sim::SimChannel> find_simchannel_ptr(raw::ChannelID_t channel) const;
  std::vector<double> sim_ides_to_xyz(std::vector<sim::IDE> const& ides) const;

  std::vector<art::Ptr<sim::SimChannel>> simchannels_; // per event, sorted by channel
  std::map<std::string, ViewOffsets> offsets_;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv_;
  mutable std::set<std::string> warned_tool_types_;
};

//-----------------------------------------------------------------------
class TriggerAnaTree : public art::EDAnalyzer {
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

  using TriggerPrimitiveWriter = VectorFieldsBuffer<TriggerPrimitiveRow>;
  using TriggerPrimitiveBacktrackingWriter = VectorFieldsBuffer<TriggerPrimitiveBacktrackingRow>;
  using TriggerPrimitiveAssociationWriter = VectorFieldsBuffer<TriggerPrimitiveAssociationRow>;
  using TriggerActivityWriter = VectorFieldsBuffer<dunedaq::trgdataformats::TriggerActivityData>;
  using TriggerActivityAssociationWriter = VectorFieldsBuffer<TriggerActivityAssociationRow>;
  using TriggerCandidateWriter = VectorFieldsBuffer<dunedaq::trgdataformats::TriggerCandidateData>;

  // Writers bound to one output tree, per tree kind
  struct TPWriters {
    TriggerPrimitiveWriter             tp;
    TriggerPrimitiveBacktrackingWriter bt;
    TriggerPrimitiveAssociationWriter  assn;
    void enable(bool backtracking, bool association) { bt.enable(backtracking); assn.enable(association); }
    void clear() { tp.clear(); bt.clear(); assn.clear(); }
    void make_branches(TTree &t) { tp.make_branches(t); bt.make_branches(t); assn.make_branches(t); }
  };
  struct TAWriters {
    TriggerActivityWriter            ta;
    TriggerActivityAssociationWriter assn;
    void enable(bool /*backtracking*/, bool association) { assn.enable(association); }
    void clear() { ta.clear(); assn.clear(); }
    void make_branches(TTree &t) { ta.make_branches(t); assn.make_branches(t); }
  };
  struct TCWriters {
    TriggerCandidateWriter tc;
    void enable(bool /*backtracking*/, bool /*association*/) {}
    void clear() { tc.clear(); }
    void make_branches(TTree &t) { tc.make_branches(t); }
  };
  template <typename Writers>
  struct TreeSet {
    TTree  *tree = nullptr;
    Writers writers;
  };

  art::ServiceHandle<art::TFileService> tfs;

  std::unordered_map<int, int> trkid_to_truth_block_id;
  std::unordered_map<int, std::string> truth_block_id_to_generator_name;
  std::map<std::string, TreeSet<TPWriters>> tp_trees; // key: raw tag
  std::map<std::string, TreeSet<TAWriters>> ta_trees;
  std::map<std::string, TreeSet<TCWriters>> tc_trees;
  std::map<int, double> track_en_sums;
  std::map<int, double> track_electron_sums;
  // map for tracking true visible energy deposited on each apa rop (for ROI studies).

  struct TPCEnergyData
  {
    double energy = 0.;
    double num_electrons = 0.;
  };
  
  std::map<ChannelInfo, TPCEnergyData> simide_tpc_energy_map;
  std::set<int> simchannel_tpcsets_; // TPCSets covered by the matched SimChannel collections (per event)
  std::map<art::ProductID, std::string> tp_tool_type_cache_; // per event


  bool dump_tp = false, dump_ta = false, dump_tc = false;
  art::InputTag tp_tag_regex, ta_tag_regex, tc_tag_regex;

  bool tp_backtracking = false;
  bool need_truth_maps = false;

  template <typename Writers>
  TreeSet<Writers> &get_or_create_tree(std::map<std::string, TreeSet<Writers>> &trees,
                                       const std::string &dir_name, const std::string &dir_title,
                                       const std::string &tag, bool backtracking, bool association);

  static ChannelInfo get_channel_info_for_channel(geo::WireReadoutGeom const *geom, raw::ChannelID_t channel);

  // Per-product fill functions called by analyze(), in this order
  void build_truth_maps(art::Event const &e);
  void fillMCTruth(art::Event const &e);
  void fillSimChannels(art::Event const &e, geo::WireReadoutGeom const *geom);
  void fillSimIDESummary();
  void fillMCParticles(art::Event const &e);
  void fillTPs(art::Event const &e, geo::WireReadoutGeom const *geom);
  void fillTAs(art::Event const &e, geo::WireReadoutGeom const *geom);
  void fillTCs(art::Event const &e);

  ChannelInfo fill_tp_row(TPWriters &tpw,
                          const dunedaq::trgdataformats::TriggerPrimitive &tp,
                          geo::WireReadoutGeom const *geom);

  // Backtracks the TP whose row is in w.tp.row (already filled by fill_tp_row):
  // checks that the TP's TPCSet is covered by the SimChannel collections, matches
  // IDEs, fills and commits w.bt. Throws if no SimChannel collection covers the TPCSet.
  void backtrack_tp(TPWriters &w, const ChannelInfo &chinfo,
                    const std::string &tool_type, const std::string &tag);

  // tpalg.tool_type of the producer of the TP collection `id` (cached per
  // event). Throws if it cannot be determined or is not a TPC TP algorithm.
  const std::string &tp_tool_type_for(art::Event const &e, art::ProductID id,
                                      const std::string &ta_tag);

  // Event meta data buffer  
  ScalarFieldsBuffer<EventMetaData> ev_sbuf;

  // visible energy for the event
  TTree *summary_tree = nullptr;
  ScalarFieldsBuffer<EventSummaryData> evsummary_buf;

  // MCTruth
  bool dump_mctruths = false;

  TTree* mctruth_tree = nullptr;
  VectorFieldsBuffer<MCTruthRow> mctruth_buffer;

  TTree* mcneutrino_tree = nullptr;
  VectorFieldsBuffer<MCNeutrinoRow> mcneutrino_buffer;

  bool dump_mcparticles = false;

  TTree* mcparticle_tree = nullptr;
  VectorFieldsBuffer<MCParticleRow> mcparticle_buffer;

  std::unique_ptr<TPBacktracker> tp_bt_;

  bool dump_simides = false;
  art::InputTag simchannel_tag;
  TTree* simide_tree = nullptr;
  VectorFieldsBuffer<SimIDERow> simide_buffer;

  TTree* simide_summary_tree = nullptr;
  ScalarFieldsBuffer<SimIDESummaryRow> simide_summary_buffer;
  VectorFieldsBuffer<SimIDETPCRow> simide_tpc_buffer;

  // JSON metadata
  nlohmann::json info_data;
  bool mctruth_map_warned_ = false;
};

} // namespace dunetrigger

#endif // DUNETRIGGER_TRIGGERANATREE_MODULE_HH
