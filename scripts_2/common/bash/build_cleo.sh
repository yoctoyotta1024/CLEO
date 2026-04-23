#!/bin/bash

### Please note: script may assume required CLEO_[XXX]
### variables have already exported (!)

set -e

COMMON_BASH_SRC="${CLEO_PATH2CLEO}/scripts_2/common/bash/src"
MACHINE_BASH_SRC="${CLEO_PATH2CLEO}/scripts_2/${CLEO_MACHINE}/bash/src"

### -------------- prepare to build CLEO --------------- ###
if [ -f "${MACHINE_BASH_SRC}/build_flags.sh" ]; then
  source ${MACHINE_BASH_SRC}/build_flags.sh
fi

if [ -f "${COMMON_BASH_SRC}/build_yac.sh" ]; then
  source ${COMMON_BASH_SRC}/build_yac.sh
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

### ensure these directories exist (it's a good idea for later use)
mkdir -p ${CLEO_PATH2BUILD}/tmp
mkdir -p ${CLEO_PATH2BUILD}/bin
mkdir -p ${CLEO_PATH2BUILD}/share
### ---------------------------------------------------- ###
