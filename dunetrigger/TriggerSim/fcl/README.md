# TriggerSim fcl files

## Layout

```
triggersim_common_cfg.fcl    detector-independent: services_triggersim, source_triggersim, outputs_triggersim
hd/
  trigger_{primitive,activity,candidate}_makers_hd_cfg.fcl   HD TP/TA/TC maker prologs
  triggersim_hd_cfg.fcl      HD producer chains (producers_/makers_/physics_triggersim_*)
  triggersim_<geom>_<tp>_<ta>_<tc>.fcl                       HD jobs (1x2x2, 1x2x6, protodunehd)
  production/tpg_dune10kt_<geom>.fcl                         HD TP-only production jobs (SimpleThreshold + AbsRunningSum)
vd/
  trigger_{primitive,activity}_makers_vd_cfg.fcl             VD TP/TA maker prologs
  triggersim_vd_cfg.fcl      VD producer chains (producers_/makers_/physics_triggersim_vd_*)
  triggersim_vd_<geom>_<tp>_<ta>.fcl                         VD jobs
  triggersim_tpg_vd_1x8x14.fcl                               old name of triggersim_vd_1x8x14_simpleThr_swift.fcl
```

## Rules

- `*_cfg.fcl` files contain only a prolog. Every other `.fcl` is a job you can run with `lar -c`.
- Includes go in one direction only:

  ```
  job -> triggersim_<det>_cfg.fcl -> triggersim_common_cfg.fcl
                                  -> trigger_*_makers_<det>_cfg.fcl
  ```

  Maker cfg files include nothing. Jobs include their detector's `triggersim_<det>_cfg.fcl` (plus
  services/tools), never maker cfg files directly. Nothing in `vd/` includes anything in `hd/`, and
  vice versa: a VD file that refers to an HD setting must fail to parse instead of silently running
  with HD values.
- Maker prologs don't know the module labels. Each `producers_triggersim_*` chain picks the labels and
  sets every consumer's `tp_tag`/`ta_tag` to match; `makers_triggersim_*` is the matching trigger path.
  Labels name the algorithm (`tpmakerTPCSimpleThreshold`, `tamakerTPCSWIFT`, ...), and they end up in
  product and TriggerAnaTree tree names.
- Geometry-specific settings (services, overrides) go only in the job files.
- fcl files install into one flat directory, so `hd/` and `vd/` are only for organisation and every
  file name must be unique on its own: detector-specific files carry `_hd`/`_vd` or a geometry tag.

## Running fewer steps

The example jobs run every step in their chain (TP -> TA -> TC for HD, TP -> TA for VD). To run fewer,
drop the unwanted producers from the job's `physics.producers` and `physics.makers`.
