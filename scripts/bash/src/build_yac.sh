#!/bin/bash

set -e
bashsrc=${CLEO_PATH2CLEO}/scripts/bash/src

### -------------------- check inputs ------------------ ###
source ${bashsrc}/check_inputs.sh
check_args_not_empty "${CLEO_PATH2CLEO}" "${CLEO_CXX_COMPILER}" "${CLEO_ENABLEYAC}"

if [ [ ${CLEO_ENABLEYAC} != "true" || ${CLEO_YACYAXTROOT} == "" ]]
then
  echo "Bad inputs, YAC must be enabled and yacyaxtroot directory must be specified for YAC build"
  exit 1
fi

if  [ "${CLEO_CXX_COMPILER}" == "intel" ]
then
  echo "YAC build currently not compatible with intel compiler" # TODO(CB): fix this incompatibility
  exit 1
fi
### ---------------------------------------------------- ###

### ------------------ choose YAC build ---------------- ###
yac_netcdf=netcdf-c/4.8.1-openmpi-4.1.2-gcc-11.2.0 # must match gcc compiler and openmpi wrapper (!)
yac_openblas=openblas@0.3.18%gcc@=11.2.0            # must match gcc compiler (!)

module load ${yac_netcdf}
spack load ${yac_openblas}
export CLEO_YAC_FLAGS="-DENABLE_YAC_COUPLING=ON -DYAXT_ROOT=${CLEO_YACYAXTROOT}/yaxt -DYAC_ROOT=${CLEO_YACYAXTROOT}/yac"
export CLEO_MODULE_PATH="${CLEO_MODULE_PATH} ${CLEO_PATH2CLEO}/libs/coupldyn_yac/cmake"
### ---------------------------------------------------- ###
