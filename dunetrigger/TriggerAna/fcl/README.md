TEMPORARY NOTES DURING DEVELOPMENT

These fcls inherit from `triggersim.fcl`, and the three modules `TPC_InfoDisplay`, `TPC_InfoComparator` and `AnaTree` are treated in a "parallel" way.
You can always create a fcl with multiple analyzers from different modules running at the same time.

There are no producers run here, therefore these will work on a file produced from a TriggerSim fcl. 
But again, one can always create a config by adding them (examples will be in `TriggerSim/fcl`).

There are only very basic fcls, check the wholw configuration by dumping before running blindly.