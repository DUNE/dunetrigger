# TriggerTools

Utilities for trigger studies that are not part of the TriggerSim/TriggerAna chains: an event filter,
a raw-waveform analyzer and a geometry dump executable.

## Layout

```
HydraFilter_module.cc     art filter: keep events with energy deposits in a chosen TPC set (APA/CRM)
RawDigitTreeDumper_module.cc     art analyzer: ADC histograms and (optionally) waveforms from RawDigits
exe/geo_config.h          shared geometry shortcuts and wire readout selection for the executables
exe/geo_inspector.cxx     standalone executable: print a detector geometry and save it as JSON
exe/channel_dumper.cxx    standalone executable: dump the channel-to-wire mapping of a TPC/plane or TPC set/ROP
fcl/
  hydra_filter.fcl                                   HydraFilter job (VD 1x8x14)
  save_vd_<geom>_3view_30deg_rawdigit_hists.fcl       RawDigitTreeDumper jobs, histograms only
  save_vd_<geom>_3view_30deg_rawdigit_waveforms.fcl   RawDigitTreeDumper jobs, histograms + waveform tree
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

## RawDigitTreeDumper (analyzer)

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
geo_inspector <geometry> [wire_readout]
```

`<geometry>` is either a geometry table from `geometry_dune.fcl` (e.g. `dunevd10kt_1x8x6_3view_30deg_geo`)
or one of the shortcuts:

| Shortcut    | Geometry                             |
|-------------|--------------------------------------|
| `vd_1x8x6`  | `dunevd10kt_1x8x6_3view_30deg_geo`   |
| `vd_1x8x14` | `dunevd10kt_1x8x14_3view_30deg_geo`  |
| `hd_1x2x6`  | `dune10kt_1x2x6_geo`                 |
| `hd_10kt`   | `dune10kt_geo`                       |
| `vd_10kt`   | `dunevd10kt_geo`                     |
| `protodunevd` | `protodunevd_v5_geo`               |

The shortcuts and the wire readout selection are shared with `channel_dumper` through `exe/geo_config.h`;
running either executable without valid arguments (`-h` for `channel_dumper`) lists them. The wire
readout table is chosen from the geometry (CRP-based `protodunevd*` and full `dunevd10kt*` need their own
`*_wire_readout` table, which also configures the `PDMapTool`) and printed as `Wire readout: <label>`;
the optional second argument overrides it.

It prints the geometry summary, then for cryostat 0 the cryostat bounds, each TPC (center, x/y/z range,
number of channels and first channel), the number of TPC sets, and each optical detector (center,
width, height, length). The cryostat, TPC and OpDet information is also written to `<geometry>.json`
in the current directory.

## channel_dumper (executable)

Dumps the channel-to-wire mapping of one TPC (or plane) or of one TPC set (or readout plane, ROP),
for checking channel numbering and for building channel maps in trigger configurations.

Run it inside the dunesw environment, e.g. `dsw -C <area> exec channel_dumper hd_1x2x6 --tpc 1 --plane 2`.

```
channel_dumper <geo> (--tpc T [--plane P] | --tpcset S [--rop R]) [--cryo C] [--format none|csv|json|both] [-o <prefix>] [-h]
```

| Option             | Description |
|--------------------|-------------|
| `<geo>`            | Geometry table from `geometry_dune.fcl` or one of the `geo_inspector` shortcuts (`vd_1x8x6`, `vd_1x8x14`, `hd_1x2x6`, `hd_10kt`) |
| `--cryo C`         | Cryostat ID (default 0) |
| `--tpc T`          | TPC mode: all planes of TPC `T`, or only `--plane P` |
| `--tpcset S`       | ROP mode: all ROPs of TPC set `S`, or only `--rop R` |
| `--wire-readout L` | Wire readout table `L` from `geometry_dune.fcl` (default: chosen from `<geo>`, see below) |
| `--table`          | Print a channel vs `z0` table to stdout after the summary (default off) |
| `--format`         | `none` (default, stdout summary only), `csv`, `json` or `both` |
| `-o`, `--output`   | Output file prefix |
| `-h`, `--help`     | Print usage and exit |

The shortcuts are the same as for `geo_inspector` (see `-h`).

The wire readout table is chosen from the geometry name, because the CRP-based geometries need a
table that also configures the `PDMapTool`: `protodunevd*` uses `protodunevd_wire_readout`,
`dunevd10kt*` without `1x8x`/`1x6x` uses `dunevd10kt_wire_readout`, `dune35t*` uses
`dune35t_wire_readout`, everything else `dune_wire_readout`. The geometry and readout labels are printed
at startup; `--wire-readout` overrides the choice.

Exactly one of `--tpc` / `--tpcset` is required; `--plane` needs `--tpc` and `--rop` needs `--tpcset`.

In TPC mode every wire of the selected plane(s) is mapped to its channel. In ROP mode the channels
`FirstChannelInROP .. +Nchannels-1` of each selected ROP are mapped back to wires. Either way the wires
are grouped by channel: in HD the induction wires wrap around the APA, so one channel reads several
wires (possibly in different TPCs). Wires are written as `cryostat:tpc:plane:wire`.

The stdout summary gives the detector name, the selection, the number of channels and wires, the
channel range, channels per view and the number of channels with more than one wire. It then prints a table with one
row per group, one group per plane in TPC mode (also with `--plane`) and one per ROP in ROP mode:

```
  plane  view  n_ch  ch_min  ch_max  pitch [cm]  dz [cm]  z0_min [cm]  z0_max [cm]
  -----  ----  ----  ------  ------  ----------  -------  -----------  -----------
      0     U   800       0     799      0.4667   0.5748        1.464      661.296
      1     V   800     800    1599      0.4667   0.5748     -430.380      228.877
      2     Z   480    2080    2559      0.4790   0.4790        0.598      230.039
```

The first column is `rop` in ROP mode, and column widths follow the content. In ROP mode the summary also
has a `TPCs in set` line (`TPCsetToTPCs`) after `Selection`, and the group table gets a last, left-aligned
`tpcs` column with the TPCs of each ROP (`ROPtoTPCs`), e.g. `0,1`. The JSON has `selection.tpcs` for the set
and `tpcs` in each ROP group.

With `--tpcset`, the summary of each TPC of the set follows the set summary (`--- TPC T of TPC set S ---`).
These are TPC-mode summaries (per plane) built from the selected channels only (so `--rop` restricts them
too), counting only the wires that sit in that TPC. They are not written to the CSV or JSON files. An `n_nonfinite` column
appears only if some group has non-finite `z0` values.

Each row gives the number of channels, their channel range, the wire pitch and the range of `z0` over all wires of
those channels, where `z0 = WireGeo::ComputeZatY0()` is the wire's z at y = 0 (an extrapolation, so it may
lie outside the detector). In TPC mode wires are grouped by their own plane, so the wrapped HD wires
that belong to a channel but sit in another TPC still count towards `z0`. Non-finite `z0` values (wires
parallel to z) are skipped in the range and counted in `n_nonfinite`. The pitch is `PlaneGeo::WirePitch()` in cm: in
TPC mode that of the selected TPC's plane, in ROP mode that of every plane the ROP's wires belong to
(including wrapped wires in other TPCs). If those differ by more than 1e-4 cm the pitch column shows
`min-max` and a `pitch by plane:` line under the row lists each `C:T:P=pitch`. The same groups are written to
`summary.groups` in the JSON (`plane` or `rop`, `view`, `n_channels`, `channel_min`, `channel_max`,
`z0_min`, `z0_max`, `n_non_finite_z0`, `pitch`; ROP mode adds `pitch_min`, `pitch_max` and
`planes: [{plane_id, pitch}]`). `pitch` is the smallest pitch of the group. `dz` is the
wire spacing along z, `pitch / |sin(ThetaZ)|` with `WireGeo::ThetaZ()` the wire angle to the z axis,
computed per wire from its own plane's pitch; it equals the pitch for wires perpendicular to z and is larger
for angled wires (about 1.23 x pitch for HD U/V, 2 x pitch for VD U/V). It is shown as `min-max` if not
uniform in the group, and written as `dz` (smallest value) in the JSON groups. Wires parallel to z
(`|sin| < 1e-6`) give `nan` and are ignored in the range. Channels are only printed to stdout with `--table`.

`--table` prints one row per selected channel, sorted by channel, with the plane (TPC mode) or ROP
(ROP mode), the number of wires, `i_z`, `ch-i_z`, `ch+i_z` and the `z0` (cm) of each wire in the channel's wire order.
`i_z = (z0 - z0_min) / dz` is the wire index along z, rounded to the nearest integer (`std::lround`) in this
table only; `z0_min` is taken from the channel's group in the summary table and `dz` is computed for the
wire itself. `ch-i_z` is `channel - i_z` and `ch+i_z` is `channel + i_z` (signed integers, from the rounded
`i_z`): `ch-i_z` is constant where channel number and z increase together, and `ch+i_z` is constant where
they run in opposite directions (the other one then changes by 2 per channel). For single-wire channels both are shown in
their columns; for multi-wire channels they show `-` and the last column lists the z0 values, then
`| i_z`, `| ch-i_z` and `| ch+i_z` with one value per wire, then the wire IDs in parentheses, so each value can be
matched to its wire. Non-finite values print as `nan`. Neighbouring wires of a plane are `dz` apart in z0,
not `pitch` apart, for angled wires; dividing by `dz` makes `i_z` step by about 1 per channel on all
planes. In ROP mode wires from several TPCs share one `z0_min`, so their `i_z` need not be a clean index.
The CSV and JSON files are not affected.

```
  channel  view  rop  n_wires  i_z  ch-i_z  ch+i_z  z0 [cm]
  -------  ----  ---  -------  ---  ------  ------  -------
      800     V    1        3    -       -       -  661.017, -200.473, 201.204 | i_z 1899, 400, 1099 | ch-i_z -1099, 400, -299 | ch+i_z 2699, 1200, 1899 (0:0:1:0, 0:1:1:400, 0:0:1:800)
     1600     Z    2        1    0    1600    1600  0.598
```

Output files are `<prefix>.csv` and/or `<prefix>.json`. The default prefix is
`<geoid>_tpc<T>[_plane<P>]` or `<geoid>_tpcset<S>[_rop<R>]`, where `<geoid>` is the resolved geometry
table name.

- CSV: header `channel,view,n_wires,wires` (ROP mode adds a trailing `rop` column), wires joined by `;`:
  `800,V,3,0:0:1:0;0:1:1:400;0:0:1:800,1`
- JSON: `{"detector", "selection": {"mode", "cryo", "tpc"/"plane" or "tpcset"/"rop"}, "summary": {...},
  "channels": [{"channel", "view", "wires": [...]}]}`; ROP mode adds `"rop"` to each channel.

Exit codes: 0 on success or `-h`, 1 for command-line errors (message plus usage on stderr), 2 if the
cryostat, TPC, plane, TPC set or ROP does not exist.

## Known issues

- Three of the `save_vd_*` files set `rawdigit_tag` instead of `rawdigits_tag`. The parameter is
  ignored and the default `tpcrawdecoder:daq.*` is used; for the 1x8x6 files this differs from the
  intended `tpcrawdecoder:daq`.
