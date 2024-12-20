#!/bin/bash
#SBATCH --job-name=build_cleo
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=940M
#SBATCH --time=00:05:00
#SBATCH --mail-user=clara.bayley@mpimet.mpg.de
#SBATCH --mail-type=FAIL
#SBATCH --account=bm1183
#SBATCH --output=./build/bin/build_cleo_out.%j.out
#SBATCH --error=./build/bin/build_cleo_err.%j.out

### Please note: script may assume required CLEO_[XXX]
### variables have already exported (!)

set -e
module purge
spack unload --all

bashsrc=${CLEO_PATH2CLEO}/scripts/bash/src

### -------------------- check inputs ------------------ ###
if [[ "${CLEO_BUILDTYPE}" == "" || "${CLEO_COMPILERNAME}" == "" || "${CLEO_ENABLEDEBUG}" == "" ||
      "${CLEO_PATH2CLEO}" == "" || "${CLEO_PATH2BUILD}" == "" || "${CLEO_ENABLEYAC}" == "" ]]
then
  echo "Bad inputs, please check all the required inputs have been specified"
  exit 1
fi

if [[ "${CLEO_PATH2CLEO}" == "${CLEO_PATH2BUILD}" ]]
then
  echo "Bad inputs, build directory cannot match the path to CLEO source"
  exit 1
fi

if [ "${CLEO_BUILDTYPE}" != "serial" ] &&
   [ "${CLEO_BUILDTYPE}" != "openmp" ] &&
   [ "${CLEO_BUILDTYPE}" != "threads" ] &&
   [ "${CLEO_BUILDTYPE}" != "cuda" ];
then
  echo "Bad inputs, build type must be 'serial', 'openmp', 'threads' or 'cuda'"
  exit 1
fi

if [ ${CLEO_ENABLEYAC} == "true" && ${CLEO_YACYAXTROOT} == "" ]
then
  echo "Bad inputs, yacyaxtroot directory must be specified if YAC is enabled"
  exit 1
fi
### ---------------------------------------------------- ###

### -------------- prepare to build CLEO --------------- ###
source ${bashsrc}/build_basic.sh

if [[ "${buildtype}" == "openmp" ]];
then
  source ${bashsrc}/build_openmp.sh
fi
elif [[ "${buildtype}" == "threads" ]];
then
  source ${bashsrc}/build_threads.sh
fi
elif [[ "${buildtype}" == "cuda" ]];
then
  source ${bashsrc}/build_openmp.sh
  source ${bashsrc}/build_cuda.sh
fi

if [ ${CLEO_ENABLEYAC} == "true" ]
then
  source ${bashsrc}/build_yac.sh
else
  export CLEO_YAC_FLAGS="-DENABLE_YAC_COUPLING=OFF"
fi
### ---------------------------------------------------- ###

### ---------------- build CLEO with cmake ------------- ###
echo "### --------------- Build Inputs -------------- ###"
echo "CLEO_BUILDTYPE = ${CLEO_BUILDTYPE}"
echo "CLEO_COMPILERNAME = ${CLEO_COMPILERTYPE}"
echo "CLEO_PATH2CLEO = ${CLEO_PATH2CLEO}"
echo "CLEO_PATH2BUILD = ${CLEO_PATH2BUILD}"

echo "CLEO_CXX_COMPILER = ${CLEO_CXX_COMPILER}"
echo "CLEO_CC_COMPILER = ${CLEO_CC_COMPILER}"
echo "CLEO_CXX_FLAGS: ${CLEO_CXX_FLAGS}"

echo "CLEO_KOKKOS_BASIC_FLAGS: ${CLEO_KOKKOS_BASIC_FLAGS}"
echo "CLEO_KOKKOS_HOST_FLAGS: ${CLEO_KOKKOS_HOST_FLAGS}"
echo "CLEO_KOKKOS_DEVICE_FLAGS: ${CLEO_KOKKOS_DEVICE_FLAGS}"

echo "CLEO_ENABLEYAC = ${CLEO_ENABLEYAC}"
echo "CLEO_YACYAXTROOT = ${CLEO_YACYAXTROOT}"
echo "CLEO_YAC_FLAGS = ${CLEO_YAC_FLAGS}"
echo "CLEO_MODULE_PATH = ${CLEO_MODULE_PATH}"
echo "### ------------------------------------------- ###"

cmake -DCMAKE_CXX_COMPILER=${CLEO_CXX_COMPILER} \
    -DCMAKE_C_COMPILER=${CLEO_CC_COMPILER} \
    -DCMAKE_CXX_FLAGS="${CLEO_CXX_FLAGS}" \
    -DCMAKE_MODULE_PATH=${CLEO_MODULE_PATH} \
    -S ${CLEO_PATH2CLEO} -B ${CLEO_PATH2BUILD} \
    ${CLEO_KOKKOS_BASIC_FLAGS} ${CLEO_KOKKOS_HOST_FLAGS} ${CLEO_KOKKOS_DEVICE_FLAGS} \
    ${CLEO_YAC_FLAGS}

# ensure these directories exist (it's a good idea for later use)
mkdir -p ${CLEO_PATH2BUILD}/tmp
mkdir -p ${CLEO_PATH2BUILD}/bin
mkdir -p ${CLEO_PATH2BUILD}/share
### ---------------------------------------------------- ###
