# SWIFT ROP TA maker: detector parameters

Source: `src/TriggerActivityROPMakerSWIFT.cpp`.

## Can be read from services

Neither parameter is set in the HD or VD fcl, so both detectors run with the HD defaults.

| Parameter | Default | Source |
|---|---|---|
| `wire_pitch` [cm] | 0.48 (HD collection) | `geo::WireReadoutGeom::Plane(...).WirePitch()`, per ROP |
| `cm_per_tick` [cm / DTS tick] | 0.016 µs × 0.16 cm/µs | `DetectorPropertiesService::DriftVelocity()` × DTS tick period (16 ns) |

- Both set the clustering coordinates (`.cpp:138`), and so what `epsilon` means.
- `channel × wire_pitch` assumes neighbouring channels are neighbouring wires on one plane.
  - That is not true for wrapped HD induction wires.
  - Fixing it needs wire positions from the geometry.

## Built on the TPC sampling period

TP times are in DTS ticks, using the hard-coded `TPAlgTPCTool::ADC_SAMPLING_RATE_IN_DTS = 32`.

- 32 = TPC tick (0.5 µs, from `DetectorClocksService` `TPCClock().TickPeriod()`) / DTS tick (16 ns).
- `window_length` = 32000, i.e. 1000 samples × 32.
- `min_samples_over_threshold`:
  - HD: 256, i.e. 8 samples × 32.
  - **VD: 10.** That is less than one sample, so the cut does nothing. It looks like a sample count was entered by mistake.

## Detector-dependent, but tuning (stay in fcl)

- `min_adc_peak`
- `inspect_energy_threshold_sadc`, `accept_energy_threshold_sadc`, `cluster_energy_cut_sadc`
- `epsilon`, `min_samples`

## Proposed fix

In `TriggerROPActivityMakerTPC::buildAlgConfigs`, fill in `wire_pitch` and `cm_per_tick` for each ROP from the services, unless the fcl sets them.

- `tpctriggeralgs` stays free of art dependencies.
- The clock and detector-properties services need event data, so the values must be set at the first event or at begin-run, not in the constructor.
