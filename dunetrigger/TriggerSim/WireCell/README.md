# TPG inside the Wire-Cell simulation graph

Runs the TPC trigger primitive generation (the `TPAlgTPC*` tools) on the raw
digitized waveforms inside the Wire-Cell Toolkit (WCT) simulation graph of the
`tpcrawdecoder` (WireCellToolkit) module. The waveforms never become
`raw::RawDigit`s and are never written to the art/ROOT file; only the TPs are.

Currently set up for the VD 1x8x14 3view 30deg geometry.

## What it does

The dune-vd WCT graph simulates each CRM in its own pipeline
(`... -> noise -> Digitizer -> signal processing`). This package adds a
`wclsTPGFilter` node to every CRM pipeline directly after the `Digitizer`, so it
sees the raw ADC (signal + noise + baseline), before any noise filtering or
signal processing:

```
DepoSetFanout -+- per CRM: [hydra switch] ductor -> reframer -> noise -> Digitizer -> wclsTPGFilter -> (SP) -+- FrameFanin -> ...
               +- ...                                                                    |                   +-
                                                                                         v
                                                                                  wclsTPGSaver --(visit)--> art::Event
```

* `wclsTPGFilter` (one per CRM) builds each channel's waveform exactly as
  `wclsFrameSaver` builds a `raw::RawDigit` and calls the configured
  `TPAlgTPCTool` on it. The tool is the same art tool that
  `TriggerPrimitiveMakerTPC` uses (loaded with `art::make_tool`), called
  through its plane-aware `process_waveform` overload, so it never touches
  art services from the TBB worker threads. Channel planes come from the WCT
  anode.
* `wclsTPGSaver` (one per graph) collects the TPs of all CRMs (thread safe),
  sorts them by `(time_start, channel, ...)` and puts them into the event as a
  single `std::vector<dunedaq::trgdataformats::TriggerPrimitive>`, product
  `tpcrawdecoder` (no instance name).

The result is identical to running `TriggerPrimitiveMakerTPC` on RawDigits
saved from the same simulation (see Validation).

## Files

| File | |
|---|---|
| `TPGFilter.{h,cxx}`, `TPGSaver.{h,cxx}` | the WCT components, library `WireCellDuneTrigger` |
| `cfg/dune-vd/wcls-sim-drift-simchannel-tpg.jsonnet` | the graph: dunereco's `dune-vd/wcls-sim-drift-simchannel-nf-sp.jsonnet` plus the switches below (changes marked `dunetrigger:`). Installed to `wire-cell-cfg/dunetrigger/dune-vd/`, which the dunetrigger UPS table prepends to `WIRECELL_PATH` |
| `fcl/wirecell_tpg_vd_cfg.fcl` | `tpcrawdecoder` module configurations |
| `fcl/detsim_wctpg_*.fcl` | detsim jobs |

## Jobs

| fcl | TPC output |
|---|---|
| `detsim_wctpg_dunevd10kt_1x8x14_3view_30deg.fcl` | TPs + SimChannels. No signal processing, no recob::Wires. Cheapest |
| `detsim_wctpg_sp_dunevd10kt_1x8x14_3view_30deg.fcl` | TPs + SimChannels + standard SP recob::Wires (`tpcrawdecoder:gauss` etc.), can go through reco |
| `detsim_wctpg_validate_dunevd10kt_1x8x14_3view_30deg.fcl` | validation, see below |

They are the standard `standard_detsim_dunevd10kt_1x8x14_3view_30deg.fcl`
(photon detector simulation etc. unchanged) with `tpcrawdecoder` replaced.
Downstream TA makers take `tp_tag: "tpcrawdecoder"`.

## Configuration

Jsonnet switches (FHiCL `wcls_main.structs`):

| | |
|---|---|
| `use_tpg` | add the TPG nodes |
| `run_sp` | keep signal processing and the `spsignals` output. If false the TPG node drops the waveforms right after use (`pass_through: false`) |
| `tpg` | `wclsTPGFilter` configuration, Jsonnet code: `{tpalg: {...}, detid: 3, default_timestamp: 0}` |
| `use_hydra` | (from dunereco) skip CRMs without energy deposits entirely |

The `wcls_main.outputers` list must match: `wclsTPGSaver:tpgsaver` always,
`wclsFrameSaver:spsignals` only with `run_sp: true`.

**The TP algorithm is configured twice** in `wirecell_tpg_vd_cfg.fcl`
(`wctpg_vd_simpleThr`): once as Jsonnet code for the WCT graph (`wct`), and
once as a normal FHiCL `tpalg` table at the top level of the `tpcrawdecoder`
module. The FHiCL table is not used by the simulation; it is there so that
`TriggerAnaTree` (which reads the producer's `tpalg` from provenance) works.
Keep the two in sync.

`detid` defaults to 3 (`kHD_TPC`), as `TriggerPrimitiveMakerTPC` hardcodes;
`kVD_BottomTPC`/`kVD_TopTPC` are 10/11.

## Caveats

* With `use_hydra: true` (the default), CRMs without energy deposits are not
  simulated at all, so they produce no noise TPs. For TP rate / fake trigger
  studies set `use_hydra: false` (every CRM simulated: much slower).
* As for RawDigit-based TPG, each event's readout window is processed
  independently (pedestal estimate restarts every event).

## Validation

`detsim_wctpg_validate_dunevd10kt_1x8x14_3view_30deg.fcl` runs, in one job,
the graph with both the TPG nodes and the per-CRM RawDigit savers
(`tpcrawdecoder:daqtpc<N>`, `process_mode: full-sim-sp`), plus
`TriggerPrimitiveMakerTPC` (`tpmakerTPCSimpleThreshold`) on those RawDigits.
`tpcrawdecoder` and `tpmakerTPCSimpleThreshold` must contain the same TPs.
