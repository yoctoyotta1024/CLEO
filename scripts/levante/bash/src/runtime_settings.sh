#!/bin/bash

set -e
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
bashsrc=${SCRIPT_DIR}

stacksize_limit=${1} # kB

### -------------------- check inputs ------------------ ###
source ${bashsrc}/check_inputs.sh
check_args_not_empty "${stacksize_limit}" "${CLEO_BUILDTYPE}"
check_args_not_empty "${CLEO_COMPILERNAME}" "${CLEO_YACYAXTROOT}"
### ---------------------------------------------------- ###

### --------------- YAC runtime settings --------------- ###
source ${bashsrc}/levante_packages.sh
if [ "${CLEO_COMPILERNAME}" == "gcc" ]
then
  module load ${levante_gcc_openmpi}
  fyamllib=${levante_gcc_fyamllib}
elif [ "${CLEO_COMPILERNAME}" == "intel" ]
then
  module load ${levante_intel_openmpi}
  fyamllib=${levante_intel_fyamllib}
else
  echo "Bad inputs, YAC only compatible with "intel" or "gcc" compiler names"
  exit 1
fi
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${fyamllib}
export PYTHONPATH=${PYTHONPATH}:${CLEO_YACYAXTROOT}/yac/python # path to YAC python bindings
### ---------------------------------------------------- ###


### --------------- set runtime optimisations----------- ###
if [ "${CLEO_BUILDTYPE}" == "cuda" ]
then
  export UCX_RNDV_SCHEME=put_zcopy                        # Preferred communication scheme with Rendezvous protocol
  export UCX_RNDV_THRESH=16384                            # Threshold when to switch transport from TCP to NVLINK [3]
  export UCX_IB_GPU_DIRECT_RDMA=yes                       # Allow remote direct memory access from/to GPU
  export UCX_TLS=cma,rc,mm,cuda_ipc,cuda_copy,gdr_copy    # Include cuda and gdr based transport layers for communication [4]
  export UCX_MEMTYPE_CACHE=n                              # Prevent misdetection of GPU memory as host memory [5]
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
