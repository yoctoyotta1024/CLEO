#!/bin/bash

### Please note: script may assume required CLEO_[XXX]
### variables have already exported (!)

set -e
[ -f /etc/profile ] && source /etc/profile

build_cleo() {
  local common_dir="${CLEO_PATH2CLEO}/scripts_2/common"
  local machine_dir="${CLEO_PATH2CLEO}/scripts_2/${CLEO_MACHINE}"

  source "${common_dir}/check_inputs.sh"
  check_machine

  ### -------------- prepare to build CLEO --------------- ###
  if [ -f "${machine_dir}/build_flags.sh" ]; then
    source "${machine_dir}/build_flags.sh"
    if declare -F configure_machine_build_flags >/dev/null; then
      configure_machine_build_flags
    fi
  fi

  if [ -f "${machine_dir}/helpers/build_yac.sh" ]; then
    source "${machine_dir}/helpers/build_yac.sh"
    if declare -F configure_machine_yac_flags >/dev/null; then
      configure_machine_yac_flags
    fi
  fi

  ### ---------------------------------------------------- ###

  ### ---------------- build CLEO with cmake ------------- ###
  echo "### --------------- Build Flags -------------- ###"

  echo "CLEO_KOKKOS_BASIC_FLAGS: ${CLEO_KOKKOS_BASIC_FLAGS}"
  echo "CLEO_KOKKOS_HOST_FLAGS: ${CLEO_KOKKOS_HOST_FLAGS}"
  echo "CLEO_KOKKOS_DEVICE_FLAGS: ${CLEO_KOKKOS_DEVICE_FLAGS}"

  echo "CLEO_BUILD_FLAGS: ${CLEO_BUILD_FLAGS}"
  echo "CLEO_YAC_FLAGS: ${CLEO_YAC_FLAGS}"
  echo "### ------------------------------------------- ###"

  cmake -DCMAKE_CXX_COMPILER=${CLEO_CXX_COMPILER} \
      -DCMAKE_C_COMPILER=${CLEO_CC_COMPILER} \
      -DCMAKE_CXX_FLAGS="${CLEO_CXX_FLAGS}" \
      -S ${CLEO_PATH2CLEO} -B ${CLEO_PATH2BUILD} \
      ${CLEO_KOKKOS_BASIC_FLAGS} ${CLEO_KOKKOS_HOST_FLAGS} ${CLEO_KOKKOS_DEVICE_FLAGS} \
      ${CLEO_BUILD_FLAGS} ${CLEO_YAC_FLAGS}
}

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  build_cleo "$@"
fi
