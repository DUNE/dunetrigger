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

#include "FieldNames.hh"
#include "MiniBackTracker.hh"
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

constexpr int INVALID_NUM = -99999;
constexpr char INVALID_STR[] = "undef";

namespace dunetrigger {

//-----------------------------------------------------------------------
struct ChannelInfo {
  unsigned int rop_id;
  int view;
  unsigned int tpcset_id;

  bool operator<(const ChannelInfo& other) const {
    if (rop_id != other.rop_id) return rop_id < other.rop_id;
    return tpcset_id < other.tpcset_id; }
};

//-----------------------------------------------------------------------
struct EventMetaData {
  int event = INVALID_NUM;
  int run = INVALID_NUM;
  int subrun = INVALID_NUM;
};

REGISTER_FIELD_NAMES(EventMetaData,
                            event,
                            run,
                            subrun);

//-----------------------------------------------------------------------
struct EventSummaryData {
  int mctruths_count = INVALID_NUM;
  int mcparticles_count = INVALID_NUM;
  int mcneutrinos_count = INVALID_NUM;
  int simides_count = INVALID_NUM;
  double tot_visible_energy_rop0 = 0.;
  double tot_visible_energy_rop1 = 0.;
  double tot_visible_energy_rop2 = 0.;
  double tot_visible_energy_rop3 = 0.;
  double tot_numelectrons_rop0 = 0.;
  double tot_numelectrons_rop1 = 0.;
  double tot_numelectrons_rop2 = 0.;
  double tot_numelectrons_rop3 = 0.;
};

REGISTER_FIELD_NAMES(EventSummaryData,
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
                            tot_numelectrons_rop3);

//-----------------------------------------------------------------------
struct MCTruthRow {
  int pdg = INVALID_NUM;
  std::string process = INVALID_STR;
  int status_code = INVALID_NUM;
  int block_id = INVALID_NUM;
  int truth_track_id = INVALID_NUM;
  std::string generator_name = INVALID_STR;
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
};

REGISTER_FIELD_NAMES(MCTruthRow,
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
                         kinetic_energy);

//-----------------------------------------------------------------------
struct MCNeutrinoRow {
  int block_id = INVALID_NUM;
  std::string generator_name = INVALID_STR;
  int nupdg = INVALID_NUM;
  int leptonpdg = INVALID_NUM;
  int ccnc = INVALID_NUM;
  int mode = INVALID_NUM;
  int interactionType = INVALID_NUM;
  int target = INVALID_NUM;
  int hitnuc = INVALID_NUM;
  int hitquark = INVALID_NUM;
  double w = 0.;
  double x = 0.;
  double y = 0.;
  double qsqr = 0.;
  double pt = 0.;
  double theta = 0.;
};

REGISTER_FIELD_NAMES(MCNeutrinoRow,
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
                         theta);

//-----------------------------------------------------------------------
struct MCParticleRow {
  int pdg = INVALID_NUM;
  std::string generator_name = INVALID_STR;
  int status_code = INVALID_NUM;
  int g4_track_id = INVALID_NUM;
  int mother = INVALID_NUM;
  int truth_block_id = INVALID_NUM;
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
  std::string process = INVALID_STR;
};

REGISTER_FIELD_NAMES(MCParticleRow,
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
                         process);

//-----------------------------------------------------------------------
struct SimIDERow {
  unsigned int channel = 0;
  int timestamp = INVALID_NUM;
  float numelectrons = 0.f;
  float energy = 0.f;
  float x = 0.f;
  float y = 0.f;
  float z = 0.f;
  int trackID = INVALID_NUM;
  float origTrackID = 0.f;
  float readout_plane_id = 0.f;
  float readout_view = 0.f;
  float detector_element = 0.f;
};

REGISTER_FIELD_NAMES(SimIDERow,
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
                         detector_element);

//-----------------------------------------------------------------------
struct SimIDESummaryRow {
  double total_visible_energy = 0.;
  double total_numelectrons = 0.;
};

REGISTER_FIELD_NAMES(SimIDESummaryRow,
                         total_visible_energy,
                         total_numelectrons);


//-----------------------------------------------------------------------
struct SimIDETPCRow {
  int readout_plane_id = INVALID_NUM;
  int detector_element = INVALID_NUM;
  double energy_per_tpc = 0.;
  double numelectrons_per_tpc = 0.;
};

REGISTER_FIELD_NAMES(SimIDETPCRow,
                         readout_plane_id,
                         detector_element,
                         energy_per_tpc,
                         numelectrons_per_tpc);

//-----------------------------------------------------------------------
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
};

REGISTER_FIELD_NAMES(TriggerPrimitiveRow,
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
                         TPCSetID);

//-----------------------------------------------------------------------
struct TriggerPrimitiveBacktrackingRow {
  int bt_primary_track_id = INVALID_NUM;
  double bt_primary_track_numelectron_frac = INVALID_NUM;
  double bt_primary_track_energy_frac = INVALID_NUM;
  double bt_edep = INVALID_NUM;
  double bt_numelectrons = INVALID_NUM;
  double bt_x = INVALID_NUM;
  double bt_y = INVALID_NUM;
  double bt_z = INVALID_NUM;
  double bt_primary_x = INVALID_NUM;
  double bt_primary_y = INVALID_NUM;
  double bt_primary_z = INVALID_NUM;
  int bt_truth_block_id = INVALID_NUM;
  std::string bt_generator_name = INVALID_STR;

  void populate_backtracking_info(const std::vector<sim::IDE> &ides,
                                  const std::unordered_map<int, int> &trkid_to_truth_block,
                                  const std::unordered_map<int, std::string> &truth_id_to_gen,
                                  const MiniBackTracker& bt
                                );
};

REGISTER_FIELD_NAMES(TriggerPrimitiveBacktrackingRow,
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
                         bt_generator_name);

//-----------------------------------------------------------------------
struct TriggerPrimitiveAssociationRow {
  int ta_number = INVALID_NUM;
};

REGISTER_FIELD_NAMES(TriggerPrimitiveAssociationRow,
                         ta_number);


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

  art::ServiceHandle<art::TFileService> tfs;
  std::map<std::string, TTree *> tree_map;

  size_t m_tc_number = 0;  // TC index bound to the "TCnumber" ROOT branch in TA-in-TC trees

  std::vector<art::Handle<std::vector<simb::MCTruth>>> mctruth_handles;
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
  std::map<int, std::shared_ptr<MiniBackTracker>> bt_map;


  bool dump_tp, dump_ta, dump_tc;
  std::string tp_tag_regex, ta_tag_regex, tc_tag_regex;

  bool tp_backtracking;
  bool need_truth_maps;

  void make_tp_tree_if_needed(std::string tag, bool assn = false);
  void make_ta_tree_if_needed(std::string tag, bool assn = false);
  void make_tc_tree_if_needed(std::string tag);

  std::vector<sim::IDE> match_simides_to_tps(const TriggerPrimitiveRow &tp, const std::string &tool_type, const MiniBackTracker& bt) const;

  ChannelInfo get_channel_info_for_channel(geo::WireReadoutGeom const *geom, int channel);

  // Per-product fill functions called by analyze(), in this order
  void build_truth_maps(art::Event const &e);
  void fillMCTruth(art::Event const &e);
  void fillSimChannels(art::Event const &e, geo::WireReadoutGeom const *geom);
  void fillSimIDESummary();
  void fillMCParticles(art::Event const &e);
  void fillTPs(art::Event const &e, geo::WireReadoutGeom const *geom);
  void fillTAs(art::Event const &e, geo::WireReadoutGeom const *geom);
  void fillTCs(art::Event const &e);

  ChannelInfo fill_tp_row(TriggerPrimitiveWriter &tpw,
                          const dunedaq::trgdataformats::TriggerPrimitive &tp,
                          geo::WireReadoutGeom const *geom);

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
  art::InputTag simchannel_tag;
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
