#!/bin/bash

### Configures CLEO with cmake.
### Requires CLEO_* variables (incl. CLEO_MACHINE) to already be exported.

set -e
[ -f /etc/profile ] && source /etc/profile

build_cleo() {
  local common_dir
  common_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
  source "${common_dir}/check_inputs.sh"
  check_machine
  local machine_dir="${common_dir}/../${CLEO_MACHINE}"

  ### -------------- prepare to build CLEO --------------- ###
  source "${machine_dir}/build_flags.sh"
  configure_machine_build_flags

  source "${machine_dir}/helpers/build_yac.sh"
  configure_machine_yac_flags
  ### ---------------------------------------------------- ###

  ### ---------------- build CLEO with cmake ------------- ###
  echo "### --------------- Build Flags -------------- ###"
  echo "CLEO_CXX_COMPILER: ${CLEO_CXX_COMPILER}"
  echo "CLEO_CC_COMPILER: ${CLEO_CC_COMPILER}"
  echo "CLEO_CXX_FLAGS: ${CLEO_CXX_FLAGS}"
  echo "CLEO_KOKKOS_BASIC_FLAGS: ${CLEO_KOKKOS_BASIC_FLAGS}"
  echo "CLEO_KOKKOS_HOST_FLAGS: ${CLEO_KOKKOS_HOST_FLAGS}"
  echo "CLEO_KOKKOS_DEVICE_FLAGS: ${CLEO_KOKKOS_DEVICE_FLAGS}"
  echo "CLEO_BUILD_FLAGS: ${CLEO_BUILD_FLAGS}"
  echo "CLEO_YAC_FLAGS: ${CLEO_YAC_FLAGS}"
  echo "CLEO_ENABLEDEBUG: ${CLEO_ENABLEDEBUG}"
  echo "### ------------------------------------------- ###"

  if [ "${CLEO_ENABLEDEBUG}" = "true" ]; then
    CLEO_CMAKE_BUILD_TYPE="Debug"
  else
    CLEO_CMAKE_BUILD_TYPE="Release"
  fi

  cmake -DCMAKE_CXX_COMPILER=${CLEO_CXX_COMPILER} \
      -DCMAKE_C_COMPILER=${CLEO_CC_COMPILER} \
      -DCMAKE_CXX_FLAGS="${CLEO_CXX_FLAGS}" \
      -DCMAKE_BUILD_TYPE=${CLEO_CMAKE_BUILD_TYPE} \
      -S ${CLEO_PATH2CLEO} -B ${CLEO_PATH2BUILD} \
      ${CLEO_KOKKOS_BASIC_FLAGS} ${CLEO_KOKKOS_HOST_FLAGS} ${CLEO_KOKKOS_DEVICE_FLAGS} \
      ${CLEO_BUILD_FLAGS} ${CLEO_YAC_FLAGS}
  ### ---------------------------------------------------- ###
}

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  build_cleo "$@"
fi
