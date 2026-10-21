# Vendored lardata code (temporary)

`ArtDataHelper/` contains copies of `lardata/ArtDataHelper/GetManyByRegexTag.h`
and `lardata/ArtDataHelper/RegexTagSelector.{h,cxx}`, which are not yet part of
a LArSoft release.

- Source: `alessandrothea/lardata`, branch `feature/getmanybyregex_explicitprocess`,
  commit `8ec849769e3d`. It replaces the earlier `RegexTagMatcher` (vendored from
  `feature/getmanyregex_bugfix`, commit `079e218a7d20`) with an `art::SelectorBase`,
  `RegexTagSelector`, passed to `getInputTags()`. Semantics differ from the old
  matcher:
  - an empty label or instance pattern matches only an empty field (use `".*"`
    for "any");
  - an empty process keeps only the newest process with a *valid* matching
    product, falling back to older processes like `e.getHandle<T>()`;
  - an explicit process is a regex and keeps the products of *every* matching
    process (e.g. `".*"` returns duplicates from all processes).

  See the doc comment in `GetManyByRegexTag.h`.
- Changes: a two-line provenance comment at the top of each file, and
  `#include "lardata/ArtDataHelper/..."` rewritten to
  `#include "dunetrigger/vendor/lardata/ArtDataHelper/..."`.
  Code and namespace (`lar::util`) are unchanged.
- `RegexTagSelector.cxx` is built as the library
  `dunetrigger_vendor_lardata_ArtDataHelper`.

## Users

| User | Include | `CMakeLists.txt` |
|---|---|---|
| `TriggerAna/TriggerAnaTree_module.cc` | `dunetrigger/vendor/lardata/ArtDataHelper/GetManyByRegexTag.h` | `TriggerAna/CMakeLists.txt`: `dunetrigger_vendor_lardata_ArtDataHelper` in `MODULE_LIBRARIES` (`lardata::ArtDataHelper` is already linked there) |
| `TriggerTools/RawDigitAna_module.cc` | same | `TriggerTools/CMakeLists.txt` |
| `TriggerTools/SimChannelTPCDump_module.cc` | same | `TriggerTools/CMakeLists.txt` |
| `TriggerSim/TriggerPrimitiveMakerTPC_module.cc` | same | `TriggerSim/CMakeLists.txt` |

`grep -rn vendor/lardata dunetrigger` lists them.

## Removal, once a LArSoft release ships these files (with the fix)

The build stops with an error as soon as the lardata in use provides
`lardata/ArtDataHelper/RegexTagSelector.h` (guard in `ArtDataHelper/CMakeLists.txt`).

**Check the fix first.** lardata `v10_01_00` already ships `GetManyByRegexTag.h`
with the older `RegexTagMatcher.{h,cxx}`, *without* the fix: its `match()` builds
the regexes from the candidate tag (`tag.instance()`, ...) instead of the
configured pattern, so a pattern such as `"tpcrawdecoder:simpleSC.*"` matches
nothing. That release does not trip the guard (no `RegexTagSelector.h`), but do
not switch to it. Remove the vendored copy only once the installed lardata
ships `RegexTagSelector`:

```
ls $LARDATA_DIR/source/lardata/ArtDataHelper/RegexTagSelector.cxx
```

(`$LARDATA_DIR` from `ups active`.) Linking both libraries would define the
same `lar::util` symbols twice.

### TriggerAnaTree (`dunetrigger/TriggerAna/`)

1. `TriggerAnaTree_module.cc`: replace
   `#include "dunetrigger/vendor/lardata/ArtDataHelper/GetManyByRegexTag.h"` with
   `#include "lardata/ArtDataHelper/GetManyByRegexTag.h"`. It is used in
   `fillSimChannels` (`simchannel_tag`), `fillTPs`, `fillTAs`, `fillTCs`
   (`tp_tag_regex`, `ta_tag_regex`, `tc_tag_regex`). No code change: the API
   (`lar::util::getManyByRegexTag<T>(e, art::InputTag)`) is the same.
2. `TriggerAna/CMakeLists.txt`: delete the `dunetrigger_vendor_lardata_ArtDataHelper`
   line from `MODULE_LIBRARIES`; `lardata::ArtDataHelper` is already there.
3. Rebuild and rerun the reference job
   (`lar -c triggerana_vd_1x8x14_simpleThr.fcl ...`); the output must be
   identical to the previous one (`compare_trees.py`, 0 differences). A matcher
   without the fix shows up immediately as
   `Found no std::vector<sim::SimChannel> collections matching simchannel_tag`.
4. Related: once larsim >= v10_27_00 is in the stack, `TPBacktracker` can switch to `cheat::BackTrackerService`; see the next section.

### Switching `TPBacktracker` to `cheat::BackTrackerService` (needs larsim >= v10_27_00)

`TPBacktracker` (in `TriggerAna/TriggerAnaTree_module.hh/.cc`) keeps its own per-event SimChannel index (one merged
vector of all collections matched by `simchannel_tag`, sorted by channel) because
`cheat::BackTracker::PrepSimChannels` in larsim <= v10_23_00 (and v11_00_00rc0) loads a
single collection with `getValidHandle(SimChannelModuleLabel)`, while the DUNE
detsim writes one `tpcrawdecoder:simpleSC<N>` collection per TPCSet. larsim
`release/v10_27_00` (commit 058e33c, "Introducing support for single TPC detsim
products using getManyByRegexTag") loads every collection matching
`SimChannelModuleLabel` as a regex and merges them, exactly as `TPBacktracker` does.

Do not enable `BackTrackerService` before that larsim is in the stack: the old
service would fail its own `getValidHandle("tpcrawdecoder:simpleSC")`.

Once it is:

1. fcl (`TriggerAna/fcl/triggerana_common.fcl` or the job fcls):
   `services.BackTrackerService: @local::dunefd_backtrackerservice` and
   `services.BackTrackerService.BackTracker.SimChannelModuleLabel: "tpcrawdecoder:simpleSC.*"`
   (the same regex as the module's `simchannel_tag`, which stays for the SimIDE dump).
2. The `TPBacktracker` class in `TriggerAna/TriggerAnaTree_module.hh/.cc`: the SimChannel code is confined to two private
   functions. Replace the `std::vector<art::Ptr<sim::SimChannel>> simchannels_`
   member with `art::ServiceHandle<cheat::BackTrackerService> bt_serv_;`
   (`#include "larsim/MCCheater/BackTrackerService.h"`), then
   - `find_simchannel_ptr(channel)` -> `return bt_serv_->provider()->FindSimChannelPtr(channel);`
   - `sim_ides_to_xyz(ides)` -> `return bt_serv_->SimIDEsToXYZ(ides);`
   - delete `prep_event`/`clear_event` and their two call sites in
     `TriggerAnaTree_module.cc` (`fillSimChannels`, `analyze`). The service
     prepares itself at the start of every event.
   `ParticleInventoryService` is already required; nothing else changes.
3. `TriggerAna/CMakeLists.txt`: add `larsim::MCCheater_BackTrackerService_service`
   to `MODULE_LIBRARIES` in `TriggerAna/CMakeLists.txt`.
4. Rebuild, rerun the reference job; `compare_trees.py` against the previous
   output must report 0 differences. Keep the TPCSet checks in the module
   (`simchannel_tpcsets_`: overlap and coverage); the service does not detect
   duplicate channels across collections.

### Other users and the directory

1. `RawDigitAna_module.cc`, `SimChannelTPCDump_module.cc` and `TriggerPrimitiveMakerTPC_module.cc`: same include
   replacement; in `TriggerTools/CMakeLists.txt` and `TriggerSim/CMakeLists.txt`
   replace `dunetrigger_vendor_lardata_ArtDataHelper` with `lardata::ArtDataHelper`
   (or drop it where that is already linked).
2. Delete `dunetrigger/vendor/` and the
   `add_subdirectory(vendor/lardata/ArtDataHelper)` line in
   `dunetrigger/CMakeLists.txt`.
