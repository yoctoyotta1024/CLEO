#!/bin/bash

### Sets the runtime environment for running CLEO on jupiter.
### Usage: configure_machine_runtime_settings [stacksize_limit (kB)]
### Requires CLEO_BUILDTYPE, CLEO_COMPILERNAME and CLEO_YACYAXTROOT to be exported.

set -e

configure_machine_runtime_settings() {
  local stacksize_limit="${1:-}"  # kB

  local machine_dir
  machine_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
  local common_dir="${machine_dir}/../common"

  ### -------------------- check inputs ------------------ ###
  source "${common_dir}/check_inputs.sh"
  check_args_not_empty "${CLEO_BUILDTYPE}" "${CLEO_COMPILERNAME}" "${CLEO_YACYAXTROOT}"
  ### ---------------------------------------------------- ###

  ### ------------ load compiler/runtime stack ------------ ###
  source "${machine_dir}/helpers/jupiter_packages.sh"
  jupiter_load_runtime_stack "${CLEO_COMPILERNAME}" "${CLEO_BUILDTYPE}"
  jupiter_load_python "${CLEO_COMPILERNAME}"
  ### ---------------------------------------------------- ###

  ### --------------- YAC runtime settings --------------- ###
  local fyamllib
  fyamllib=$(jupiter_fyamllib_for_compiler "${CLEO_COMPILERNAME}")
  export LD_LIBRARY_PATH="${fyamllib}:${LD_LIBRARY_PATH}"
  export PYTHONPATH="${PYTHONPATH}:${CLEO_PATH2CLEO}/examples/exampleplotting/plotcleo:${CLEO_YACYAXTROOT}/yac/python"
  ### ---------------------------------------------------- ###

  ### ------------ communication runtime (MPI) ------------ ###
  if [[ "${CLEO_BUILDTYPE}" == "cuda" ]]; then
    if [[ -z "${CLEO_CUDA_ROOT}" ]]; then
      echo "Error: CLEO_CUDA_ROOT is not set for cuda runtime."
      exit 1
    fi
    export LD_LIBRARY_PATH="${CLEO_CUDA_ROOT}/lib64:${LD_LIBRARY_PATH}"
    export UCX_RNDV_SCHEME=put_zcopy                        # Preferred communication scheme with Rendezvous protocol
    export UCX_RNDV_THRESH=16384                            # Threshold when to switch transport from TCP to NVLINK
    export UCX_IB_GPU_DIRECT_RDMA=yes                       # Allow remote direct memory access from/to GPU
    export UCX_TLS=cma,rc,mm,cuda_ipc,cuda_copy,gdr_copy    # Include cuda and gdr based transport layers
    export UCX_MEMTYPE_CACHE=n                              # Prevent misdetection of GPU memory as host memory
  else
    export UCX_TLS="shm,rc_mlx5,rc_x,self" # for jobs using LESS than 150 nodes
  fi
  export OMPI_MCA_osc="ucx"
  export OMPI_MCA_pml="ucx"
  export OMPI_MCA_btl="self"
  export OMPI_MCA_pml_ucx_opal_mem_hooks=1
  export UCX_HANDLE_ERRORS="bt"
  export OMPI_MCA_io="romio321"           # basic optimisation of I/O
  ### ---------------------------------------------------- ###

  ### ------------------ threading ----------------------- ###
  export OMP_PROC_BIND=spread
  export OMP_PLACES=threads
  ### ---------------------------------------------------- ###

  ### ------------------ process limits ------------------ ###
  if [[ -n "${stacksize_limit}" ]]; then
    ulimit -s "${stacksize_limit}"
  fi
  ulimit -c 0

  # Prevent glibc from automatically trimming the heap.
  export MALLOC_TRIM_THRESHOLD_="-1"
  ### ---------------------------------------------------- ###
}

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  configure_machine_runtime_settings "$@"
fi
