# dunetrigger

Tools for emulating dunedaq trigger flow and analyzing trigger information in the DUNE experiment within the LArSoft framework. It's noteworthy that the LArSoft trigger emulation software is still under development and the commands presented here are subject to changes.

## Apptainer container

Currently, the apptainer container is used to run the LArSoft trigger emulation software. In this [link](https://wiki.dunescience.org/wiki/SL7_to_Alma9_conversion#Running_SL7_inside_a_container), one can find the instructions to set up the apptainer environment, which is a containerized environment for running `dunesw` tools on SL7 instead of using AL9, still under development.

To sum up, the following command can be used to run the apptainer container:

```bash
/cvmfs/oasis.opensciencegrid.org/mis/apptainer/current/bin/apptainer shell --shell=/bin/bash -B /cvmfs,/exp,/nashome,/pnfs/dune,/opt,/run/user,/etc/hostname,/etc/hosts,/etc/krb5.conf --ipc --pid /cvmfs/singularity.opensciencegrid.org/fermilab/fnal-dev-sl7:latest
```

## Create and set up the development environment

The bash snippet below creates and sets up the development environment for the LArSoft trigger emulation software used in the DUNE experiment. To do so, one just needs to copy the code into an executable bash script and run it.

It is worth noting that the username must be adapted to the user running the script. Regarding the software version, it is set to `v09_89_01d01`, but it can be changed to any other future version available in the DUNE software environment. The same applies to the qualifiers, which are set to `e26:prof` in this case.

```bash
# This script sets up the development environment for the LArSoft trigger emulation software used in the DUNE experiment.
# It creates a working directory, sets up the required software packages, and clones the necessary repositories.
# Finally, it builds and installs the software.

VERSION=v09_89_01d01  # Main version of the software to be used
QUALS=e26:prof  # Qualifiers for the software packages
DIRECTORY=trigger_sim_dev  # Name of the directory where the software will be installed
export WORKDIR=/exp/dune/app/users/$USER/  # Path to the working directory

# Create the working directory if it doesn't exist
if [ ! -d "$WORKDIR" ]; then
  mkdir -p "$WORKDIR"
fi

# Set up the DUNE software environment
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
export UPS_OVERRIDE="-H Linux64bit+3.10-2.17"
setup dunesw -v ${VERSION} -q ${QUALS}

cd ${WORKDIR}
mkdir -p ${DIRECTORY}
cd ${DIRECTORY}

# Create a new development area using mrb
mrb newDev -v ${VERSION} -q ${QUALS}
source ${WORKDIR}/${DIRECTORY}/localProducts*/setup

# Clone the dunetrigger repository
mrb g https://github.com/wesketchum/dunetrigger.git@develop

# Clone the duneprototypes repository
mrb g https://github.com/wesketchum/duneprototypes.git

cd ${MRB_SOURCE}/duneprototypes/
git checkout feature/wketchum_TriggerDecodingUpdates

# Clone the dunecore repository
mrb g https://github.com/wesketchum/dunecore.git

cd ${MRB_SOURCE}/dunecore/
git checkout feature/wketchum_TriggerDecodingUpdates

cd ${MRB_BUILDDIR}

# Set up the environment for building and installing the software
mrbsetenv
mrb i -j8
mrbslp
```

## Source the local products after creating the development area

Once the development area was created, the environment can be set up by running the following commands within a bash script:

```bash
# This script sets up the environment for the DUNE software and runs the necessary commands to source the local products after creating the development area.

VERSION=v09_89_01d01  # Version of the software to be used
QUALS=e26:prof  # Qualifiers for the software packages

# Source the setup script for the DUNE software
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh

# Set the UPS_OVERRIDE environment variable to specify the desired platform
export UPS_OVERRIDE="-H Linux64bit+3.10-2.17"

# Setup the specific version of the DUNE software
setup dunesw ${VERSION} -q ${QUALS}

# Source the setup script for the local products associated to the development area
source localProducts_*/setup

# Set up the MRB source local products
mrbslp
```

As mentioned in the previous [section](#create-and-set-up-the-development-environment), both the software version and qualifiers can be changed to any other available in the DUNE software environment.

## Generate trigger information from real raw data

Generating trigger information from real data is essential to validate the correct performance of trigger algorithms already available in the DAQ framework. The goal is to compare the trigger information generated in LArSoft with the trigger information obtained from the DAQ framework and stored in the *Trigger Records*.

The following subsections are devoted to generate trigger information (TPs, TAs and TCs) from raw waveforms from data taking using the LArSoft trigger emulation software.

### How to decode a raw trigger data file

The following command shows how to decode a raw trigger data file using the LArSoft trigger emulation software.

```bash
lar -c run_pdhd_tpc_decoder.fcl -n ${N_EVENTS} -s ${RAW_FILE_PATH} -o ${DECODED_FILE_PATH} -T ${DECODED_HIST_PATH}
```

Here is a brief explanation of the different flags used in the command:

- `-n ${N_EVENTS}`: specifies the number of events to be processed.
- `-s ${RAW_FILE_PATH}`: specifies the path to the raw trigger data file, an HDF5 file with trigger objects. E.g., `np04hd_raw_run026305_0033_dataflow0_datawriter_0_20240520T133910.hdf5`.
- `-o ${DECODED_FILE_PATH}`: specifies the path for the output decoded art ROOT file.
- `-T ${DECODED_HIST_PATH}`: specifies the path for the output histogram ROOT file.

The configuration file `run_pdhd_tpc_decoder.fcl` can be found in the `/srcs/dunetrigger/example` directory of the dunetrigger repository.

In case one desires to use a raw data file for tests, one can find an example raw data file in the `/pnfs/dune/persistent/users/hamza/trigger_sim_testing` directory.

For the example raw data file, assuming that the working directory is `/srcs/dunetrigger/example`, the command would be:

```bash
lar -c run_pdhd_tpc_decoder.fcl -n -1 -s /pnfs/dune/persistent/users/hamza/trigger_sim_testing/np04hd_raw_run026305_0033_dataflow0_datawriter_0_20240520T133910.hdf5
```

The `-1` flag indicates that all events in the file will be processed. Since `-o` and `-T` flags were not utilized, the output files will be saved in the working directory, having names defined, by default, in the fhicl file.

### How to run the LArSoft trigger emulation

The following command shows how to run the LArSoft trigger emulation software for the decoded file, which contains the raw waveforms needed for TPG finding.

```bash
lar -c run_tpalg_taalg_tcalg.fcl -n ${N_EVENTS} -s ${DECODED_FILE_PATH} -o ${TRIGGER_FILE_PATH} -T ${TRIGGER_HIST_PATH}
```

In this command, `${N_EVENTS}` represents the number of events to be processed. `${DECODED_FILE_PATH}` corresponds to the path of the decoded art ROOT file. `${TRIGGER_FILE_PATH}` denotes the path for the output art ROOT data file. Lastly, `${TRIGGER_HIST_PATH}` indicates the path for the output histogram ROOT file. The configuration file `run_tpalg_taalg_tcalg.fcl` can be located in the `/srcs/dunetrigger/example` directory of the dunetrigger repository.

One has to modify the `run_tpalg_taalg_tcalg.fcl` file to set the desired trigger algorithms and/or parameters for the trigger emulation. By default, the fhicl file is set to use the `TXAlgTPCExample` algorithms.

```c++
producers:
 {
     tpmakerTPC:
     {
         module_type: TriggerPrimitiveMakerTPC
         rawdigit_tag: "tpcrawdecoder:daq"
         tpalg: {
             tool_type: TPAlgTPCExample
             verbosity: 1
         }
         verbosity: 1
     }

     tamakerTPC:
     {
         module_type: TriggerActivityMakerTPC
         tp_tag: "tpmakerTPC"
         taalg: {
             tool_type: TAAlgTPCExample
             multiplicity: 100
             verbosity: 1
         }
         verbosity: 1
     }

     tcmakerTPC:
     {
         module_type: TriggerCandidateMakerTPC
         ta_tag: "tamakerTPC"
         tcalg: {
             tool_type: TCAlgTPCExample
             multiplicity: 1
             verbosity: 1
         }
         verbosity: 1
     }

 }
 ```

At present, there is a more sophisticated algorithm for trigger primitive generation available, called `TPAlgTPCSimpleThreshold`. It emulates the trigger primitive generator used for data taking during NP-04 tests and runs in 2024. This one should be implemented in the fhicl file as follows:

```c++
tpmakerTPC:
     {
         module_type: TriggerPrimitiveMakerTPC
         rawdigit_tag: "tpcrawdecoder:daq"
         tpalg: {
             tool_type: TPAlgTPCSimpleThreshold
             threshold_tpg_plane0: 500
             threshold_tpg_plane1: 500
             threshold_tpg_plane2: 500
             verbosity: 0
         }
         verbosity: 1
     }
```

One may notice that there is no actual *Trigger Activity/Candidate Maker* included in the LArSoft trigger emulation software. The `TAAlgTPCExample` and `TCAlgTPCExample` placeholders are used to exemplify the implementation of these algorithms. The user can implement their own algorithms for these modules.

### How to compare the trigger data

To check the reproducibility of the trigger data, by comparing both offline and online trigger information, one can utilize the following instruction:

```shell
lar -c run_triggerTPCInfoComparator.fcl -n ${N_EVENTS} -s ${TRIGGER_FILE_PATH} -o ${TRIGGER_COMPARATOR_FILE_PATH} -T ${TRIGGER_COMPARATOR_HIST_PATH}
```

In this instruction, `${N_EVENTS}` specifies the number of events to process. `${TRIGGER_FILE_PATH}` denotes the location of the trigger data file, while `${TRIGGER_COMPARATOR_FILE_PATH}` points to the output art ROOT data file location. `${TRIGGER_COMPARATOR_HIST_PATH}` indicates the location for the output histogram ROOT file. The `run_triggerTPCInfoComparator.fcl` configuration file can be found in the `/srcs/dunetrigger/example` directory within the dunetrigger repository.

## Generate trigger information from simulated data with LArSoft

Using the art ROOT file generated at the `detsim` level, one can generate trigger information by means of the LArSoft trigger emulation software.

### LArSoft trigger emulation from simulated data

To perform the LArSoft trigger emulation, just execute the `run_tpalg_taalg_tcalg.fcl` file as in the previous [section](#how-to-run-the-larsoft-trigger-emulation) for the decoded art ROOT file. Keep in mind that the input file, `${DETSIM_FILE_PATH}`,  is the art ROOT file generated at the `detsim` level.

```bash
lar -c run_tpalg_taalg_tcalg.fcl -n ${N_EVENTS} -s ${DETSIM_FILE_PATH} -o ${TRIGGER_FILE_PATH} -T ${TRIGGER_HIST_PATH}
```

### How to analyze the trigger data

To analyze the trigger data, either simulated or raw-decoded, one can use the following command:

```bash
lar -c run_offlineTriggerTPCInfoDisplay.fcl -n ${N_EVENTS} -s ${TRIGGER_FILE_PATH} -o ${ANALYSIS_FILE_PATH} -T ${ANALYSIS_HIST_PATH}
```

In this command, `${N_EVENTS}` is the number of events to be processed, `${TRIGGER_FILE_PATH}` is the path to the trigger data file, `${ANALYSIS_FILE_PATH}` is the path for the output art ROOT data file, and `${ANALYSIS_HIST_PATH}` is the path for the output histogram ROOT file. The configuration file `run_offlineTriggerTPCInfoDisplay.fcl` can be found in the `/srcs/dunetrigger/example` directory of the dunetrigger repository.