#ifndef DUNETRIGGER_TRIGGERANAROWS_HH
#define DUNETRIGGER_TRIGGERANAROWS_HH

#include "detdataformats/trigger/TriggerActivityData.hpp"
#include "detdataformats/trigger/TriggerCandidateData.hpp"
#include "detdataformats/trigger/TriggerPrimitive.hpp"

#include "FieldNames.hh"

#include <cstdint>
#include <string>

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
struct TriggerActivityAssociationRow {
  int tc_number = INVALID_NUM;
};

REGISTER_FIELD_NAMES(TriggerActivityAssociationRow,
                         tc_number);

// dunedaq data structs used directly as buffer rows
REGISTER_FIELD_NAMES(dunedaq::trgdataformats::TriggerActivityData,
                         version, time_start, time_end, time_peak, time_activity,
                         channel_start, channel_end, channel_peak,
                         adc_integral, adc_peak, detid, type, algorithm);

REGISTER_FIELD_NAMES(dunedaq::trgdataformats::TriggerCandidateData,
                         version, time_start, time_end, time_candidate, detid, type, algorithm);
} // namespace dunetrigger

#endif // DUNETRIGGER_TRIGGERANAROWS_HH
