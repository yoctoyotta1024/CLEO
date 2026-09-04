#!/bin/bash

set -e

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)

COMMON_BASH_SRC="${SCRIPT_DIR}/../common"
jupiter_HELPERS_DIR="${SCRIPT_DIR}/helpers"

configure_machine_runtime_settings() {

    local stacksize_limit="${1}"

    source "${COMMON_BASH_SRC}/check_inputs.sh"
    source "${jupiter_HELPERS_DIR}/jupiter_packages.sh"

    check_args_not_empty \
        "${stacksize_limit}" \
        "${CLEO_BUILDTYPE}" \
        "${CLEO_COMPILERNAME}" \
        "${CLEO_YACYAXTROOT}"

    # Load compiler/runtime and Python environment.
    jupiter_load_runtime_stack "${CLEO_COMPILERNAME}" "${CLEO_BUILDTYPE}"
    jupiter_load_python "${CLEO_COMPILERNAME}"

    # Runtime library paths.
    local fyamllib
    fyamllib=$(jupiter_fyamllib_for_compiler "${CLEO_COMPILERNAME}")

    export LD_LIBRARY_PATH="${fyamllib}:${LD_LIBRARY_PATH}"

    export PYTHONPATH="${PYTHONPATH}:${CLEO_PATH2CLEO}/examples/exampleplotting/plotcleo:${CLEO_YACYAXTROOT}/yac/python"

    # Configure communication runtime.
    if [[ "${CLEO_BUILDTYPE}" == "cuda" ]]; then

        if [[ -z "${CLEO_CUDA_ROOT}" ]]; then
            echo "Error: CLEO_CUDA_ROOT is not set for CUDA runtime."
            exit 1
        fi

        export LD_LIBRARY_PATH="${CLEO_CUDA_ROOT}/lib64:${LD_LIBRARY_PATH}"

        # GPU-aware UCX configuration.
        export UCX_RNDV_SCHEME=put_zcopy
        export UCX_RNDV_THRESH=16384
        export UCX_IB_GPU_DIRECT_RDMA=yes
        export UCX_TLS=cma,rc,mm,cuda_ipc,cuda_copy,gdr_copy
        export UCX_MEMTYPE_CACHE=n

    else

        # CPU UCX configuration.
        export UCX_TLS="shm,rc_mlx5,rc_x,self"

    fi

    # Use UCX for MPI communication.
    export OMPI_MCA_osc="ucx"
    export OMPI_MCA_pml="ucx"
    export OMPI_MCA_btl="self"
    export OMPI_MCA_pml_ucx_opal_mem_hooks=1

    export UCX_HANDLE_ERRORS="bt"
    export OMPI_MCA_io="romio321"

    # OpenMP thread placement.
    export OMP_PROC_BIND=spread
    export OMP_PLACES=threads

    # Prevent glibc from automatically trimming the heap.
    export MALLOC_TRIM_THRESHOLD_="-1"

    # Process limits.
    ulimit -s "${stacksize_limit}"
    ulimit -c 0
}

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    configure_machine_runtime_settings "$@"
fi
