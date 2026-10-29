#!/bin/bash

set -e

configure_openmp_build() {
  local common_dir
  common_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)

  ### -------------------- check inputs ------------------ ###
  source "${common_dir}/check_inputs.sh"
  check_args_not_empty "${CLEO_BUILDTYPE}"

  if [[ "${CLEO_BUILDTYPE}" != "openmp" ]]; then
    echo "Bad inputs, build type for enabling openmp on host must be 'openmp'"
    exit 1
  fi
  ### ---------------------------------------------------- ###

  ### ------- choose host parallelism kokkos flags ------- ###
  export CLEO_KOKKOS_HOST_FLAGS="${CLEO_KOKKOS_HOST_FLAGS} -DKokkos_ENABLE_OPENMP=ON"
  ### ---------------------------------------------------- ###
}

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  configure_openmp_build "$@"
fi
