#include "TPBacktracker.hh"

#include "cetlib_except/exception.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "dunetrigger/TriggerSim/TPAlgTools/TPAlgTPCTool.hh"

#include <algorithm>
#include <iterator>
#include <utility>

int dunetrigger::ViewOffsets::for_view(geo::View_t view) const {
  switch (view) {
    case geo::kU:
      return u;
    case geo::kV:
      return v;
    case geo::kW:
      return w;
    default:
      return 0;
  }
}

dunetrigger::TPBacktracker::TPBacktracker(std::map<std::string, ViewOffsets> offsets)
    : offsets_(std::move(offsets)) {}

std::vector<sim::IDE> dunetrigger::TPBacktracker::match_ides(const TriggerPrimitiveRow& tp,
                                                             const std::string& tool_type,
                                                             const MiniBackTracker& bt) const {
  auto it = offsets_.find(tool_type);
  if (it == offsets_.end() && warned_tool_types_.insert(tool_type).second) {
    mf::LogWarning("TPBacktracker") << "No offsets found for tool type " << tool_type << ", using 0,0,0";
  }
  const ViewOffsets offsets = it != offsets_.end() ? it->second : ViewOffsets{};
  int offset = offsets.for_view(static_cast<geo::View_t>(tp.readout_view));
  int sample_start = tp.time_start / TPAlgTPCTool::ADC_SAMPLING_RATE_IN_DTS;
  int sample_end = sample_start + tp.samples_over_threshold;
  sample_start += offset;
  sample_end += offset;
  sample_start = std::max(0, sample_start);
  sample_end = std::max(0, sample_end);

  if (sample_start > sample_end) {
    throw cet::exception("TPBacktracker") << "Invalid sample range";
  }

  art::Ptr<sim::SimChannel> sim_channel = bt.findSimChannelPtr(tp.channel);
  // No SimChannel for this channel (e.g. a noise-only TP): nothing to match.
  if (sim_channel.isNull()) {
    return {};
  }
  std::vector<sim::IDE> matched_ides = sim_channel->TrackIDsAndEnergies(sample_start, sample_end);
  return matched_ides;
}

void dunetrigger::TPBacktracker::fill_row(TriggerPrimitiveBacktrackingRow& row,
                                          const std::vector<sim::IDE>& ides,
                                          const std::unordered_map<int, int>& trkid_to_truth_block,
                                          const std::unordered_map<int, std::string>& truth_id_to_gen,
                                          const MiniBackTracker& bt) const {
  row.bt_primary_track_id = INVALID_NUM;
  row.bt_primary_track_numelectron_frac = INVALID_NUM;
  row.bt_primary_track_energy_frac = INVALID_NUM;
  row.bt_edep = 0.;
  row.bt_numelectrons = 0.;
  row.bt_x = INVALID_NUM;
  row.bt_y = INVALID_NUM;
  row.bt_z = INVALID_NUM;
  row.bt_primary_x = INVALID_NUM;
  row.bt_primary_y = INVALID_NUM;
  row.bt_primary_z = INVALID_NUM;
  row.bt_truth_block_id = INVALID_NUM;
  row.bt_generator_name.clear();

  if (ides.empty()) {
    return;
  }

  // Resolve each IDE (with a non-zero trackID) to its MC track once.
  std::vector<std::pair<const sim::IDE*, int>> resolved;
  resolved.reserve(ides.size());
  size_t n_unresolved = 0;
  for (const sim::IDE& ide : ides) {
    if (ide.trackID == 0) {
      continue;
    }
    const simb::MCParticle* part = pi_serv_->TrackIdToParticle_P(ide.trackID);
    if (!part) {
      ++n_unresolved;
      continue;
    }
    resolved.emplace_back(&ide, part->TrackId());
  }
  if (n_unresolved > 0) {
    mf::LogDebug("TPBacktracker") << n_unresolved << " IDEs skipped: trackID not resolved to an MCParticle";
  }

  if (resolved.empty()) {
    mf::LogDebug("TPBacktracker") << "Empty IDEs set!";
    return;
  }

  std::map<int, double> track_numelectrons;
  std::map<int, double> track_energies;

  for (const auto& [ide, mc_track_id] : resolved) {
    track_numelectrons[mc_track_id] += ide->numElectrons;
    track_energies[mc_track_id] += ide->energy;
    row.bt_numelectrons += ide->numElectrons;
    row.bt_edep += ide->energy;
  }

  row.bt_primary_track_id =
      std::max_element(track_numelectrons.begin(), track_numelectrons.end(), [](const auto& a, const auto& b) {
        return a.second < b.second;
      })->first;

  std::vector<sim::IDE> primary_ides;
  for (const auto& [ide, mc_track_id] : resolved) {
    if (mc_track_id == row.bt_primary_track_id) {
      primary_ides.push_back(*ide);
    }
  }

  row.bt_primary_track_numelectron_frac = track_numelectrons[row.bt_primary_track_id] / row.bt_numelectrons;
  row.bt_primary_track_energy_frac = track_energies[row.bt_primary_track_id] / row.bt_edep;

  std::vector<double> bt_position = bt.simIDEsToXYZ(ides);
  std::vector<double> primary_bt_position = bt.simIDEsToXYZ(primary_ides);

  row.bt_x = bt_position[0];
  row.bt_y = bt_position[1];
  row.bt_z = bt_position[2];
  row.bt_primary_x = primary_bt_position[0];
  row.bt_primary_y = primary_bt_position[1];
  row.bt_primary_z = primary_bt_position[2];

  auto tb_it = trkid_to_truth_block.find(row.bt_primary_track_id);
  if (tb_it == trkid_to_truth_block.end()) {
    mf::LogDebug("TPBacktracker") << "No truth block for primary track " << row.bt_primary_track_id;
    return;
  }
  row.bt_truth_block_id = tb_it->second;
  auto gen_it = truth_id_to_gen.find(row.bt_truth_block_id);
  if (gen_it == truth_id_to_gen.end()) {
    mf::LogDebug("TPBacktracker") << "No generator name for truth block " << row.bt_truth_block_id;
    return;
  }
  row.bt_generator_name = gen_it->second;
}
