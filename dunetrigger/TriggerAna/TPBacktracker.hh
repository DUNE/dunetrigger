#ifndef DUNETRIGGER_TRIGGERANA_TPBACKTRACKER_HH
#define DUNETRIGGER_TRIGGERANA_TPBACKTRACKER_HH

#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

#include "dunetrigger/TriggerAna/TriggerAnaRows.hh"

#include <array>
#include <map>
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

// Matches TPs to SimChannel IDEs and to MC truth. SimChannel lookup goes
// through cheat::BackTrackerService, which merges every SimChannel collection
// matching its BackTracker.SimChannelModuleLabel. That label must equal the
// analyzer's simchannel_tag; check_simchannels() verifies it every event.
class TPBacktracker {
public:
  explicit TPBacktracker(std::map<std::string, ViewOffsets> offsets);

  // Throws if the number of SimChannels held by cheat::BackTrackerService
  // differs from n_expected (the count in the collections matching
  // simchannel_tag), i.e. if services.BackTrackerService.BackTracker.
  // SimChannelModuleLabel does not match simchannel_tag.
  void check_simchannels(std::size_t n_expected, std::string const& simchannel_tag) const;

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
  // The only two places that look up SimChannels / IDE positions: thin wrappers
  // around cheat::BackTrackerService (provider()->FindSimChannelPtr / SimIDEsToXYZ).
  art::Ptr<sim::SimChannel> find_simchannel_ptr(raw::ChannelID_t channel) const;
  std::vector<double> sim_ides_to_xyz(std::vector<sim::IDE> const& ides) const;

  std::map<std::string, ViewOffsets> offsets_;
  art::ServiceHandle<cheat::BackTrackerService const> bt_serv_;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv_;
  mutable std::set<std::string> warned_tool_types_;
};

} // namespace dunetrigger

#endif // DUNETRIGGER_TRIGGERANA_TPBACKTRACKER_HH
