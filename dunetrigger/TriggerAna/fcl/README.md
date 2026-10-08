# TriggerAna fcl files

## Layout

```
triggerana_common.fcl        detector-independent: services/source/outputs_triggerana_tree and the
                             generic triggerAnaTree_dumpAll analyzer (includes triggersim_common_cfg.fcl)
hd/
  triggerana_tree.fcl        HD physics chains: TriggerSim producers + TriggerAnaTree
  triggerana_tpc_infodisplay.fcl, triggerana_tpc_infocomparator.fcl
                             HD TriggerTPCInfoDisplay / TriggerTPCInfoComparator analyzers and chains
  triggerana_tree_<geom>_*.fcl, triggerana_tpc_info*_protodunehd_*.fcl, ta_dump_dune10kt_1x2x6.fcl
                             HD jobs
  production/                dune10kt 1x2x2/1x2x6 jobs for the trigger group production (2025)
vd/
  triggerana_vd_<geom>_simpleThr.fcl   VD TriggerAnaTree jobs (include triggerana_common.fcl only)
```

The TriggerSim-only production jobs (`tpg_dune10kt_*.fcl`, TP generation with both the SimpleThreshold
and AbsRunningSum algorithms) are in `TriggerSim/fcl/hd/production/`.

## Rules

The same rules as in `TriggerSim/fcl/README.md` apply. In particular, nothing in `vd/` includes
anything in `hd/`: VD jobs include `triggerana_common.fcl`, never `triggerana_tree.fcl`.

## Notes

- The jobs either run the TriggerSim producers in the same job (`physics_triggerana_*` chains in
  `triggerana_tree.fcl` and the info display/comparator files) or read a file made by a TriggerSim job
  (analyzers only, e.g. the VD jobs and `production/triggerana_dune10kt_*.fcl`).
- The `TriggerAnaTree`, `TriggerTPCInfoDisplay` and `TriggerTPCInfoComparator` analyzers are
  independent: a job can run any combination of them.
- TriggerAnaTree names its TP/TA/TC trees after the producer labels, so a renamed label renames the
  trees. The VD jobs select the TP trees with `tp_tag_regex: ".*SimpleThreshold.*"`.
- These are basic examples: dump the configuration (`fhicl-dump`) before running one blindly.
