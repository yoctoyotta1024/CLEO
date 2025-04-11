#!/bin/bash

set -e
source /etc/profile
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
bashsrc=${SCRIPT_DIR}
cleo_yac_module_path="${CLEO_PATH2CLEO}/libs/coupldyn_yac/cmake"

### -------------------- check inputs ------------------ ###
source ${bashsrc}/check_inputs.sh
check_args_not_empty "${CLEO_PATH2CLEO}" "${CLEO_COMPILERNAME}" "${CLEO_CXX_COMPILER}" "${CLEO_ENABLEYAC}"

if [[ "${CLEO_ENABLEYAC}" != "true" || "${CLEO_YACYAXTROOT}" == "" ]]
then
  echo "Bad inputs, YAC must be enabled and yacyaxtroot directory must be specified for YAC build"
  exit 1
fi

if  [[ "${CLEO_COMPILERNAME}" == "gcc" &&
       "${CLEO_CXX_COMPILER}" != "/sw/spack-levante/openmpi-4.1.2-mnmady/bin/mpic++" ]]
then
  echo "YAC currently requires gcc/11.2.0-gcc-11.2.0 compilers"
  exit 1
elif  [ "${CLEO_COMPILERNAME}" == "intel" ]
then
  echo "YAC build currently not compatible with intel compiler" # TODO(CB): fix this incompatibility
  exit 1
fi


### ---------------------------------------------------- ###

### ------------------ choose YAC build ---------------- ###
source ${bashsrc}/levante_packages.sh
module load ${levante_gcc_netcdf_yac}
spack load ${levante_gcc_openblas_yac}
export CLEO_YAC_FLAGS="-DCLEO_YAC_MODULE_PATH="${cleo_yac_module_path}" \
  -DCLEO_YAXT_ROOT=${CLEO_YACYAXTROOT}/yaxt \
  -DCLEO_YAC_ROOT=${CLEO_YACYAXTROOT}/yac"
### ---------------------------------------------------- ###
