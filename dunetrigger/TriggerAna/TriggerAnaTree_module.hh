#ifndef DUNETRIGGER_TRIGGERANATREE_MODULE_HH
#define DUNETRIGGER_TRIGGERANATREE_MODULE_HH

#include "detdataformats/trigger/TriggerPrimitive.hpp"
#include "lardataobj/Simulation/SimChannel.h"

#include <string>
#include <unordered_map>
#include <vector>

constexpr int INVALID = -99999;

namespace dunetrigger {

struct ChannelInfo {
  unsigned int rop_id;
  int view;
  unsigned int tpcset_id;
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

} // namespace dunetrigger

#endif // DUNETRIGGER_TRIGGERANATREE_MODULE_HH
