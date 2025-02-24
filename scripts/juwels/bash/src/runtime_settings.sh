#!/bin/bash

set -e
bashsrc=${CLEO_PATH2CLEO}/scripts/juwels/bash/src

stacksize_limit=${1} # kB

### -------------------- check inputs ------------------ ###
source ${bashsrc}/check_inputs.sh
check_args_not_empty "${stacksize_limit}" "${CLEO_BUILDTYPE}" "${CLEO_ENABLEYAC}"
check_args_not_empty "${CLEO_COMPILERNAME}" "${CLEO_ENABLE_MPTRAC}"
### ---------------------------------------------------- ###

### --------------- YAC runtime settings --------------- ###
if [ "${CLEO_ENABLEYAC}" == "true" ]
then
  echo "Bad inputs, YAC build enabled but building CLEO with YAC on JUWELS is not currently supported"
  exit 1
fi
### ---------------------------------------------------- ###

### ----------- MPTRAC compile-time settings ----------- ###
if [ "${CLEO_ENABLE_MPTRAC}" == "true" ]
then
  if [ "${CLEO_COMPILERNAME}" != "gcc" ]
  then
    echo "Bad inputs, MPTRAC can only supported with gcc compilation on JUWELS"
    exit 1
  fi
  module load ${juwels_gcc_gsl}
  module load ${juwels_gcc_netcdf}
fi
### ---------------------------------------------------- ###

### --------------- set runtime optimisations----------- ###
if [ "${CLEO_BUILDTYPE}" == "cuda" ]
then
  echo "Bad inputs, CUDA build enabled but building CLEO with CUDA on JUWELS is not currently supported"
  exit 1
fi

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export SRUN_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}
export OMP_PROC_BIND=spread # (!) will be overriden by KMP_AFFINITY
export OMP_PLACES=threads # (!) will be overriden by KMP_AFFINITY

export MALLOC_TRIM_THRESHOLD_="-1"

ulimit -s ${stacksize_limit}
ulimit -c 0
### ---------------------------------------------------- ###
