#ifndef DUNETRIGGER_TRIGGERANA_TPBACKTRACKER_HH
#define DUNETRIGGER_TRIGGERANA_TPBACKTRACKER_HH

#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "lardataobj/Simulation/SimChannel.h"
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

} // namespace dunetrigger

#endif // DUNETRIGGER_TRIGGERANA_TPBACKTRACKER_HH
