#!/bin/bash

set -e
bashsrc=${CLEO_PATH2CLEO}/scripts/juwels/bash/src

stacksize_limit=${1} # kB

### -------------------- check inputs ------------------ ###
source ${bashsrc}/check_inputs.sh
check_args_not_empty "${stacksize_limit}" "${CLEO_BUILDTYPE}" "${CLEO_ENABLEYAC}"
### ---------------------------------------------------- ###

### --------------- YAC runtime settings --------------- ###
if [ "${CLEO_ENABLEYAC}" == "true" ]
then
  check_args_not_empty "${CLEO_YACYAXTROOT}"
  source ${bashsrc}/levante_packages.sh

  spack load ${levante_gcc_python_yac}
  spack load ${levante_gcc_cython_yac}
  spack load ${levante_gcc_mpi4py_yac}

  export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${levante_gcc_fyamllib}
  export PYTHONPATH=${PYTHONPATH}:${CLEO_YACYAXTROOT}/yac/python # path to YAC python bindings
fi
### ---------------------------------------------------- ###


### --------------- set runtime optimisations----------- ###
if [ "${CLEO_BUILDTYPE}" == "cuda" ]
then
  echo "Bad inputs, CUDA build enabled but building CLEO with CUDA on JUWELS is not currently supported"
  exit 1
fi

export OMPI_MCA_osc="ucx"
export OMPI_MCA_pml="ucx"
export OMPI_MCA_btl="self"
export UCX_HANDLE_ERRORS="bt"
export OMPI_MCA_pml_ucx_opal_mem_hooks=1
export OMPI_MCA_io="romio321"          # basic optimisation of I/O
export UCX_TLS="shm,rc_mlx5,rc_x,self" # for jobs using LESS than 150 nodes

export OMP_PROC_BIND=spread # (!) will be overriden by KMP_AFFINITY
export OMP_PLACES=threads # (!) will be overriden by KMP_AFFINITY
export KMP_AFFINITY="granularity=fine,scatter" # (similar to OMP_PROC_BIND=spread)
export KMP_LIBRARY="turnaround"

export MALLOC_TRIM_THRESHOLD_="-1"

ulimit -s ${stacksize_limit}
ulimit -c 0
### ---------------------------------------------------- ###
