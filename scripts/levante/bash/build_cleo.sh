#!/bin/bash
#SBATCH --job-name=build_cleo
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --gpus=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=940M
#SBATCH --time=00:05:00
#SBATCH --account=bm1183
#SBATCH --output=./build/bin/build_cleo_out.%j.out
#SBATCH --error=./build/bin/build_cleo_err.%j.out

### Please note: script may assume required CLEO_[XXX]
### variables have already exported (!)

set -e
source /etc/profile
module purge
spack unload --all

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
bashsrc=${SCRIPT_DIR}/src

### -------------------- check inputs ------------------ ###
source ${bashsrc}/check_inputs.sh
check_args_not_empty "${CLEO_BUILDTYPE}" "${CLEO_COMPILERNAME}" "${CLEO_ENABLEDEBUG}" \
  "${CLEO_PATH2CLEO}" "${CLEO_PATH2BUILD}" "${CLEO_YACYAXTROOT}"

check_source_and_build_paths
check_buildtype
check_compilername
check_yac
### ---------------------------------------------------- ###

### -------------- prepare to build CLEO --------------- ###
source ${bashsrc}/build_basic.sh
source ${bashsrc}/build_yac.sh

if [[ "${CLEO_BUILDTYPE}" == "openmp" ]];
then
  source ${bashsrc}/build_openmp.sh
elif [[ "${CLEO_BUILDTYPE}" == "threads" ]];
then
  source ${bashsrc}/build_threads.sh
elif [[ "${CLEO_BUILDTYPE}" == "cuda" ]];
then
  source ${bashsrc}/build_openmp.sh
  source ${bashsrc}/build_cuda.sh
fi
### ---------------------------------------------------- ###

### ---------------- build CLEO with cmake ------------- ###
echo "### --------------- Build Inputs -------------- ###"
echo "CLEO_BUILDTYPE: ${CLEO_BUILDTYPE}"
echo "CLEO_COMPILERNAME: ${CLEO_COMPILERNAME}"
echo "CLEO_PATH2CLEO: ${CLEO_PATH2CLEO}"
echo "CLEO_PATH2BUILD: ${CLEO_PATH2BUILD}"

echo "CLEO_CXX_COMPILER: ${CLEO_CXX_COMPILER}"
echo "CLEO_CC_COMPILER: ${CLEO_CC_COMPILER}"
echo "CLEO_CXX_FLAGS: ${CLEO_CXX_FLAGS}"

echo "CLEO_KOKKOS_BASIC_FLAGS: ${CLEO_KOKKOS_BASIC_FLAGS}"
echo "CLEO_KOKKOS_HOST_FLAGS: ${CLEO_KOKKOS_HOST_FLAGS}"
echo "CLEO_KOKKOS_DEVICE_FLAGS: ${CLEO_KOKKOS_DEVICE_FLAGS}"

echo "CLEO_BUILD_FLAGS: ${CLEO_BUILD_FLAGS}"
echo "CLEO_YACYAXTROOT: ${CLEO_YACYAXTROOT}"
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
