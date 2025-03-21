# dunetrigger

Tools for emulating dunedaq trigger flow and analyzing trigger information in the DUNE experiment within the LArSoft framework.

## Apptainer container

Currently, the apptainer container is used to run the LArSoft trigger emulation software. In this [link](https://wiki.dunescience.org/wiki/SL7_to_Alma9_conversion#Running_SL7_inside_a_container), one can find the instructions to set up the apptainer environment, which is a containerized environment for running `dunesw` tools on SL7 instead of using AL9, still under development.

To sum up, the following command can be used to run the apptainer container:

```bash
/cvmfs/oasis.opensciencegrid.org/mis/apptainer/current/bin/apptainer shell --shell=/bin/bash -B /cvmfs,/exp,/nashome,/pnfs/dune,/opt,/run/user,/etc/hostname,/etc/hosts,/etc/krb5.conf --ipc --pid /cvmfs/singularity.opensciencegrid.org/fermilab/fnal-dev-sl7:latest
```

## Set up dunetrigger

After an initial development phase, this is now a standard repository under DUNE.
Once the desired dunesw version is set up, the modules and fcls of dunetrigger can be used the standard way with the `lar` command.

Anybody is invited to contribute to the development, by installing a local version of `dunetrigger` and developing in a branch.
Instructions on how to set up a local install can be found in any larsoft tutorial; the most recent and complete one was held in 2025 at CERN, to access the [full event](https://indico.cern.ch/event/1461779/) ask [Dominic](mailto:d.brailsford@lancaster.ac.uk) for the password.

## Convert real raw data to art ROOT

Generating trigger information from real data is essential to validate the correct performance of trigger algorithms already available in the DAQ framework. The goal is to compare the trigger information generated in LArSoft with the trigger information obtained from the DAQ framework and stored in the *Trigger Records*.

The following command shows how to decode a raw trigger data file using the LArSoft trigger emulation software.

```bash
lar -c run_pdhd_tpc_decoder.fcl -n ${N_EVENTS} -s ${RAW_FILE_PATH} -o ${DECODED_FILE_PATH} -T ${DECODED_HIST_PATH}
```

Here is a brief explanation of the different flags used in the command:

- `${N_EVENTS}`: specifies the number of events to be processed.
- `${RAW_FILE_PATH}`: specifies the path to the raw trigger data file, an HDF5 file with trigger objects. E.g., `np04hd_raw_run026305_0033_dataflow0_datawriter_0_20240520T133910.hdf5`.
- `${DECODED_FILE_PATH}`: specifies the path for the output decoded art ROOT file.
- `${DECODED_HIST_PATH}`: specifies the path for the output histogram ROOT file.

In case one desires to use a raw data file for tests, one can find an example raw data file in the `/pnfs/dune/persistent/users/hamza/trigger_sim_testing` directory.

## Run the LArSoft trigger emulation

The trigger emulation can be run using the fcl files under `dunetrigger/TriggerSim/fcl`. 
- `triggersim_makers.fcl` is a PROLOG-only configuraiton containing some defaults for the different makers. 
More example configurations can be added.
- `triggersim.fcl` is another PROLOG-only configuration that contains some example blocks that are imported to all other fcls. 
- `triggersim_*_simpleThr_simpleWin_simpleWin.fcl`, where instead of * you will find different geometries, are some configurations to run just the producers that create the TPs, TAs, TCs. 
The `simpleThr_simpleWin_simpleWin` in the names refer to the algorithms used for the three stages.

You can create your own fcl inheriting the default blocks and overwriting what you like.
If you think that a configuration can be of common use, commit it to a branch and submit a Merge Request. 
Note that for it to be useful, it needs to respect the current include logic and naming scheme, look at the examples.

To run, you need to have an input file coming either from the decoding of raw data (see previous section) or from `detsim`. 
The output is a file equivalent to a `reco1` file coming from the standard simulation chain.
The command will be, for example:

```bash
lar -c triggersim_protodunehd_simpleThr_simpleWin_simpleWin.fcl -n ${N_EVENTS} -s ${INPUT_FILE_PATH} -o ${RECO1_FILE_PATH} -T ${TRIGGERSIM_HIST_PATH}
```
Running the trigger emulation on pure simulation files coming frmo `detsim` is one of the main motives behind the creation of this package; the goal is to validate samples and reconstruction algorithms against the capabilities of the DAQ, to make sure that analyses are not relying on unrealistic assumptions.
In other words, we want to make sure that we trigger on what one wants to use for their selection, to prevent surprises once operations start.

## Compare the trigger data

To compare offline and online trigger information, use the `dunetrigger/TriggerAna/fcl/triggerana_tpc_infocomparator_<geometry>_<tpalg>_<taalg>_<tcalg>.fcl` configuration files. 
This analysis helps identify any differences between the trigger data generated by the trigger emulation software and the trigger information saved in the *Trigger Records*. 
It provides histograms that highlight the similarities and discrepancies between offline and online trigger information. 
This analysis validates the correctness and consistency of the trigger algorithms in LArSoft, ensuring accurate event selection for simulated events.

An example:

```bash
lar -c triggerana_tpc_infocomparator_protodunehd_simpleThr_simpleWin_simpleWin.fcl -n ${N_EVENTS} -s ${TRIGGER_FILE_PATH} -o ${TRIGGER_COMPARATOR_FILE_PATH} -T ${TRIGGER_COMPARATOR_HIST_PATH}
```

- `${N_EVENTS}` specifies the number of events to process. 
- `${TRIGGER_FILE_PATH}` is the full path of the converted raw data file.
- `${TRIGGER_COMPARATOR_FILE_PATH}` points to the output art ROOT data file location. 
- `${TRIGGER_COMPARATOR_HIST_PATH}` indicates the location for the output histogram ROOT file.

## Analyze the trigger data

To analyze trigger data from simulated detector signals, use the `dunetrigger/TriggerAna/fcl/triggerana_tpc_infodisplay_<geometry>_<tpalg>_<taalg>_<tcalg>.fcl` configuration files. 
This analysis provides insights into the accuracy and consistency of trigger information generated by the trigger emulation software. 
It is suitable not only to test new algorithms, but also to fine-tune algorithm parameters, allowing for proper event selections.
It includes histograms highlighting trigger primitives (TPs), trigger activities (TAs), and trigger candidates (TCs). 
By examining these outputs, one can understand the performance and behaviour of trigger algorithms and ensure the reliability of the trigger emulation system.

As an example, one can use the following command:

```bash
lar -c triggerana_tpc_infodisplay_protodunehd_simpleThr_simpleWin_simpleWin.fcl -n ${N_EVENTS} -s ${TRIGGER_FILE_PATH} -o ${TRIGGER_DISPLAY_FILE_PATH} -T ${TRIGGER_DISPLAY_HIST_PATH}
```

- `${N_EVENTS}` specifies the number of events to process. 
- `${TRIGGER_FILE_PATH}` is the full path of the converted raw data file.
- `${TRIGGER_DISPLAY_FILE_PATH}` points to the output art ROOT data file location. 
- `${TRIGGER_DISPLAY_HIST_PATH}` indicates the location for the output histogram ROOT file.

## Copyright and Licensing
Copyright © 2024 FERMI NATIONAL ACCELERATOR LABORATORY for the benefit of the DUNE Collaboration.

This repository, and all software contained within, is licensed under
the Apache License, Version 2.0 (the "License"); you may not use this
file except in compliance with the License. You may obtain a copy of
the License at

    http://www.apache.org/licenses/LICENSE-2.0

Copyright is granted to FERMI NATIONAL ACCELERATOR LABORATORY on behalf
of the Deep Underground Neutrino Experiment (DUNE). Unless required by
applicable law or agreed to in writing, software distributed under the
License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for
the specific language governing permissions and limitations under the
License.

