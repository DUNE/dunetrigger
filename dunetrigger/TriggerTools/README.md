# TriggerTools

Utilities for trigger studies that are not part of the TriggerSim/TriggerAna chains: an event filter,
a raw-waveform analyzer and a geometry dump executable.

## Layout

```
HydraFilter_module.cc     art filter: keep events with energy deposits in a chosen TPC set (APA/CRM)
RawDigitAna_module.cc     art analyzer: ADC histograms and (optionally) waveforms from RawDigits
exe/geo_inspector.cxx     standalone executable: print a detector geometry and save it as JSON
fcl/
  hydra_filter.fcl                                   HydraFilter job (VD 1x8x14)
  save_vd_<geom>_3view_30deg_rawdigi_hists.fcl       RawDigitAna jobs, histograms only
  save_vd_<geom>_3view_30deg_rawdigi_waveforms.fcl   RawDigitAna jobs, histograms + waveform tree
```

## HydraFilter (filter)

Selects events with simulated activity in one TPC set (an APA for HD, a CRM for VD), e.g. to keep
only the events relevant for a single readout unit before running the slower TPC simulation steps.

An event passes if at least one `sim::SimEnergyDeposit` with `NumElectrons() > 0` has its start or end
point inside the active bounding box of any TPC belonging to TPC set `tpcset_id` of cryostat 0.

| Parameter   | Type   | Description                     |
|-------------|--------|---------------------------------|
| `tpcset_id` | uint16 | TPC set (APA/CRM) to select on  |

- The input label is hard-coded to `IonAndScint`; events without that collection are rejected.
- Per-TPC deposit counts are printed for every event.
- `fcl/hydra_filter.fcl` runs it on the VD 1x8x14 geometry with `tpcset_id: 66` and writes the
  selected events (all products kept) to `hydra_filtered_%ifb.root`.

## RawDigitAna (analyzer)

Studies raw TPC waveforms (`raw::RawDigit`) and their ADC distributions, separating signal from noise
using the truth in `sim::SimChannel`.

For each event:

1. Channels with at least one IDE in any SimChannel collection matching `simchannels_tag` are flagged
   as signal channels.
2. The ADCs of every RawDigit collection matching `rawdigits_tag` are histogrammed per plane.
   The `Noise` histograms are filled only from channels without IDEs.
3. If `save_waveforms` is set, the waveforms are stored in `rawdigis_tree`.

| Parameter                      | Type     | Default                      | Description |
|--------------------------------|----------|------------------------------|-------------|
| `rawdigits_tag`                | InputTag | `tpcrawdecoder:daq.*`        | RawDigit collections (regex allowed) |
| `simchannels_tag`              | InputTag | `tpcrawdecoder:simpleSC.*`   | SimChannel collections (regex allowed) |
| `save_waveforms`               | bool     | `true`                       | Fill `rawdigis_tree` with waveforms |
| `save_channels_with_ides_only` | bool     | `false`                      | Store only signal-channel waveforms and add the `active_channels` branch |

The job throws if no collection matches either tag.

Output (TFileService):

- `ADCsByPlane`, `ADCsNoiseByPlane`: TH2I, ADC (0–16384) vs plane index (0 = U, 1 = V, 2 = X).
- `ADCsPlane{U,V,X}`, `ADCsNoisePlane{U,V,X}`: TH1I, ADC distribution per plane.
- `rawdigis_tree`: `event`, `run`, `subrun`, one branch per channel of cryostat 0 (branch name = channel
  number, `std::vector<short>` of ADCs) and, with `save_channels_with_ides_only`, `active_channels`.
  Branches of channels that were not saved are left empty for that event. The tree is filled only
  when `save_waveforms` is true.

The `fcl/save_vd_*` files run it on the VD 1x8x6 and 1x8x14 3-view 30° geometries: the `_hists` jobs
set `save_waveforms: false` and write `trigger_digits_hists_%ifb.root`, the `_waveforms` jobs keep the
default and write `trigger_digits_waves_%ifb.root`. Waveform output is large (one branch per channel):
restrict it to a few events or use `save_channels_with_ides_only`.

## geo_inspector (executable)

Loads a DUNE geometry outside of an art job and dumps its layout, for checking TPC/channel/optical
detector numbering when writing trigger configurations or analysis code.

```
geo_inspector <geometry>
```

`<geometry>` is either a geometry table from `geometry_dune.fcl` (e.g. `dunevd10kt_1x8x6_3view_30deg_geo`)
or one of the shortcuts:

| Shortcut    | Geometry                             |
|-------------|--------------------------------------|
| `vd_1x8x6`  | `dunevd10kt_1x8x6_3view_30deg_geo`   |
| `vd_1x8x14` | `dunevd10kt_1x8x14_3view_30deg_geo`  |
| `hd_1x2x6`  | `dune10kt_1x2x6_geo`                 |
| `hd_10kt`   | `dune10kt_geo`                       |

It prints the geometry summary, then for cryostat 0 the cryostat bounds, each TPC (center, x/y/z range,
number of channels and first channel), the number of TPC sets, and each optical detector (center,
width, height, length). The cryostat, TPC and OpDet information is also written to `<geometry>.json`
in the current directory.

## Known issues

- Three of the `save_vd_*` files set `rawdigit_tag` instead of `rawdigits_tag`. The parameter is
  ignored and the default `tpcrawdecoder:daq.*` is used; for the 1x8x6 files this differs from the
  intended `tpcrawdecoder:daq`.
