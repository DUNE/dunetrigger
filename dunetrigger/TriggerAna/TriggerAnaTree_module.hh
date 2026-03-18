#ifndef DUNETRIGGER_TRIGGERANATREE_MODULE_HH
#define DUNETRIGGER_TRIGGERANATREE_MODULE_HH

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "detdataformats/trigger/TriggerActivityData.hpp"
#include "detdataformats/trigger/TriggerCandidateData.hpp"
#include "detdataformats/trigger/TriggerPrimitive.hpp"
#include "fhiclcpp/ParameterSet.h"
#include "larcore/Geometry/WireReadout.h"
#include "lardataobj/Simulation/SimChannel.h"

#include "ScalarFieldsBuffer.hh"
#include "VectorFieldsBuffer.hh"

#include <TTree.h>

#include <nlohmann/json.hpp>

#include <array>
#include <map>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

constexpr int INVALID = -99999;

namespace dunetrigger {

struct ChannelInfo {
  unsigned int rop_id;
  int view;
  unsigned int tpcset_id;

  bool operator<(const ChannelInfo& other) const {
    if (rop_id != other.rop_id) return rop_id < other.rop_id;
    return tpcset_id < other.tpcset_id; }
};

struct EventMetaData {
  int event = -1;
  int run = -1;
  int subrun = -1;
};

struct EventSummaryData {
  int mctruths_count = -1;
  int mcparticles_count = -1;
  int mcneutrinos_count = -1;
  int simides_count = -1;
  double tot_visible_energy_rop0 = 0.;
  double tot_visible_energy_rop1 = 0.;
  double tot_visible_energy_rop2 = 0.;
  double tot_visible_energy_rop3 = 0.;
  double tot_numelectrons_rop0 = 0.;
  double tot_numelectrons_rop1 = 0.;
  double tot_numelectrons_rop2 = 0.;
  double tot_numelectrons_rop3 = 0.;
};

struct MCTruthRow {
  int pdg = -1;
  std::string process = "undef";
  int status_code = -1;
  int block_id = -1;
  int truth_track_id = -1;
  std::string generator_name = "undef";
  double x = 0.;
  double y = 0.;
  double z = 0.;
  double t = 0.;
  double px = 0.;
  double py = 0.;
  double pz = 0.;
  double p = 0.;
  double energy = 0.;
  double kinetic_energy = 0.;

  MCTruthRow() = default;
};

struct MCNeutrinoRow {
  int block_id = -1;
  std::string generator_name = "undef";
  int nupdg = -1;
  int leptonpdg = -1;
  int ccnc = -1;
  int mode = -1;
  int interactionType = -1;
  int target = -1;
  int hitnuc = -1;
  int hitquark = -1;
  double w = 0.;
  double x = 0.;
  double y = 0.;
  double qsqr = 0.;
  double pt = 0.;
  double theta = 0.;

  MCNeutrinoRow() = default;
};

struct MCParticleRow {
  int pdg = -1;
  std::string generator_name = "undef";
  int status_code = -1;
  int g4_track_id = -1;
  int mother = -1;
  int truth_block_id = -1;
  double x = 0.;
  double y = 0.;
  double z = 0.;
  double t = 0.;
  double end_x = 0.;
  double end_y = 0.;
  double end_z = 0.;
  double end_t = 0.;
  double px = 0.;
  double py = 0.;
  double pz = 0.;
  double energy = 0.;
  double kinetic_energy = 0.;
  double edep = 0.;
  double numelectrons = 0.;
  double shower_edep = 0.;
  double shower_numelectrons = 0.;
  std::string process = "undef";

  MCParticleRow() = default;
};

struct SimIDERow {
  unsigned int channel = 0;
  int timestamp = -1;
  float numelectrons = 0.f;
  float energy = 0.f;
  float x = 0.f;
  float y = 0.f;
  float z = 0.f;
  int trackID = -1;
  float origTrackID = 0.f;
  float readout_plane_id = 0.f;
  float readout_view = 0.f;
  float detector_element = 0.f;

  SimIDERow() = default;
};

struct SimIDESummaryRow {
  double total_visible_energy = 0.;
  double total_numelectrons = 0.;

  SimIDESummaryRow() = default;
};


struct SimIDETPCRow {
  int readout_plane_id = -1;
  int detector_element = -1;
  double energy_per_tpc = 0.;
  double numelectrons_per_tpc = 0.;

  SimIDETPCRow() = default;
};

struct TriggerPrimitiveRow {
  uint8_t version = 0;
  uint8_t flag = 0;
  uint8_t detid = 0;
  uint32_t channel = 0;
  uint16_t samples_over_threshold = 0;
  uint64_t time_start = 0;
  uint16_t samples_to_peak = 0;
  uint32_t adc_integral = 0;
  uint16_t adc_peak = 0;
  unsigned int readout_plane_id = 0;
  int readout_view = 0;
  unsigned int TPCSetID = 0;

  void from_tp(const dunedaq::trgdataformats::TriggerPrimitive &tp);

  TriggerPrimitiveRow() = default;
};

struct TriggerPrimitiveBacktrackingRow {
  int bt_primary_track_id = INVALID;
  double bt_primary_track_numelectron_frac = INVALID;
  double bt_primary_track_energy_frac = INVALID;
  double bt_edep = 0.;
  double bt_numelectrons = 0.;
  double bt_x = INVALID;
  double bt_y = INVALID;
  double bt_z = INVALID;
  double bt_primary_x = INVALID;
  double bt_primary_y = INVALID;
  double bt_primary_z = INVALID;
  int bt_truth_block_id = INVALID;
  std::string bt_generator_name = "undef";

  void populate_backtracking_info(const std::vector<sim::IDE> &ides,
                                  const std::unordered_map<int, int> &trkid_to_truth_block,
                                  const std::unordered_map<int, std::string> &truth_id_to_gen);

  TriggerPrimitiveBacktrackingRow() = default;
};

struct TriggerPrimitiveAssociationRow {
  int ta_number = -1;

  TriggerPrimitiveAssociationRow() = default;
};

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

  art::ServiceHandle<art::TFileService> tfs;
  std::map<std::string, TTree *> tree_map;

  size_t m_tc_number = 0;  // TC index bound to the "TCnumber" ROOT branch in TA-in-TC trees

  std::unordered_map<int, int> trkId_to_truthBlockId;
  std::unordered_map<int, std::string> truthBlockId_to_generator_name;
  std::map<std::string, std::tuple<TriggerPrimitiveWriter, TriggerPrimitiveBacktrackingWriter, TriggerPrimitiveAssociationWriter>> tp_writers;

  std::map<std::string, dunedaq::trgdataformats::TriggerActivityData> ta_bufs;
  std::map<std::string, dunedaq::trgdataformats::TriggerCandidateData> tc_bufs;
  std::map<int, double> track_en_sums;
  std::map<int, double> track_electron_sums;
  // map for tracking true visible energy deposited on each apa rop (for ROI studies).

  struct TPCEnergyData
  {
    double energy;
    double num_electrons;
  };
  
  std::map<ChannelInfo, TPCEnergyData> simide_tpc_energy_map;

  bool dump_tp, dump_ta, dump_tc;
  std::string tp_tag_regex, ta_tag_regex, tc_tag_regex;

  bool tp_backtracking;

  void make_tp_tree_if_needed(std::string tag, bool assn = false);
  void make_ta_tree_if_needed(std::string tag, bool assn = false);
  void make_tc_tree_if_needed(std::string tag);

  std::vector<sim::IDE> match_simides_to_tps(const TriggerPrimitiveRow &tp, const std::string &tool_type) const;

  ChannelInfo get_channel_info_for_channel(geo::WireReadoutGeom const *geom, int channel);

  // Event meta data buffer  
  ScalarFieldsBuffer<EventMetaData> ev_sbuf;

  // visible energy for the event
  TTree *summary_tree;
  ScalarFieldsBuffer<EventSummaryData> evsummary_buf;

  // MCTruth
  bool dump_mctruths;

  TTree* mctruth_tree;
  VectorFieldsBuffer<MCTruthRow> mctruth_buffer;

  TTree* mcneutrino_tree;
  VectorFieldsBuffer<MCNeutrinoRow> mcneutrino_buffer;

  bool dump_mcparticles;

  TTree* mcparticle_tree;
  VectorFieldsBuffer<MCParticleRow> mcparticle_buffer;

  std::map<std::string, std::array<int, 3>> bt_view_offsets;

  bool dump_simides;
  std::string simchannel_tag;
  TTree* simide_tree;
  VectorFieldsBuffer<SimIDERow> simide_buffer;

  TTree* simide_summary_tree;
  ScalarFieldsBuffer<SimIDESummaryRow> simide_summary_buffer;
  VectorFieldsBuffer<SimIDETPCRow> simide_tpc_buffer;

  // JSON metadata
  nlohmann::json info_data;
  bool first_event_flag;
};

} // namespace dunetrigger

#endif // DUNETRIGGER_TRIGGERANATREE_MODULE_HH
