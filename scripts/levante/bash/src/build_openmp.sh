#!/bin/bash

set -e
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
bashsrc=${SCRIPT_DIR}

### -------------------- check inputs ------------------ ###
source ${bashsrc}/check_inputs.sh
check_args_not_empty "${CLEO_BUILDTYPE}"

if [[ "${CLEO_BUILDTYPE}" != "openmp" && "${CLEO_BUILDTYPE}" != "cuda" ]];
then
  echo "Bad inputs, build type for enabling openmp on host must be 'openmp' or 'cuda'"
  exit 1
fi
### --------------------------------------------------- ###

### ------- choose host parallelism kokkos flags ------- ###
export CLEO_KOKKOS_HOST_FLAGS="${CLEO_KOKKOS_HOST_FLAGS} -DKokkos_ENABLE_OPENMP=ON"
### ---------------------------------------------------- ###
