#ifndef DUNETRIGGER_TPBACKTRACKER_HH
#define DUNETRIGGER_TPBACKTRACKER_HH

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

#include "MiniBackTracker.hh"
#include "TriggerAnaRows.hh"

#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

namespace dunetrigger {

// Sample-window offsets applied before matching a TP to the SimChannel IDEs,
// one per readout view (U, V, W).
//
// NOTE (future improvement, not a bug): the right offsets depend on the TP
// algorithm *and its configuration* (thresholds, running-sum window, ...), i.e.
// on the TP collection, not only on the tool type. At the moment one set of
// offsets is configured per tool_type (fcl `bt_window_offsets`) and applied to
// every TP collection produced with that tool, whatever its parameters. Keying
// the offsets by TP collection, or deriving them from the algorithm parameters
// recorded in the producer's provenance, would remove this approximation.
struct ViewOffsets {
  int u = 0, v = 0, w = 0;
  int for_view(geo::View_t view) const;   // kU->u, kV->v, kW->w, else 0
};

class TPBacktracker {
public:
  explicit TPBacktracker(std::map<std::string, ViewOffsets> offsets);

  // Returns IDEs on the TP's channel inside the TP's time window shifted by the
  // tool-specific view offset. Empty if no SimChannel exists for the channel.
  std::vector<sim::IDE> match_ides(const TriggerPrimitiveRow& tp,
                                   const std::string& tool_type,
                                   const MiniBackTracker& bt) const;

  // Fills `row` from `ides`. Leaves INVALID_NUM defaults if `ides` is empty or
  // no IDE could be resolved to an MCParticle.
  void fill_row(TriggerPrimitiveBacktrackingRow& row,
                const std::vector<sim::IDE>& ides,
                const std::unordered_map<int,int>& trkid_to_truth_block,
                const std::unordered_map<int,std::string>& truth_id_to_gen,
                const MiniBackTracker& bt) const;

  const std::map<std::string, ViewOffsets>& offsets() const { return offsets_; }

private:
  std::map<std::string, ViewOffsets> offsets_;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv_;
  mutable std::set<std::string> warned_tool_types_;
};

} // namespace dunetrigger

#endif // DUNETRIGGER_TPBACKTRACKER_HH
