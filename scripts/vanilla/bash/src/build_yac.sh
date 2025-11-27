#!/bin/bash

set -e
source /etc/profile
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
bashsrc=${SCRIPT_DIR}
cleo_yac_module_path="${CLEO_PATH2CLEO}/libs/coupldyn_yac/cmake"

### -------------------- check inputs ------------------ ###
source ${bashsrc}/check_inputs.sh
check_args_not_empty "${CLEO_PATH2CLEO}" "${CLEO_COMPILERNAME}" "${CLEO_CXX_COMPILER}" "${CLEO_YACYAXTROOT}"

if  [[ "${CLEO_COMPILERNAME}" == "gcc" &&
       "${CLEO_CXX_COMPILER}" != "/sw/spack-levante/openmpi-4.1.2-mnmady/bin/mpic++" ]]
then
  echo "YAC currently requires gcc/11.2.0-gcc-11.2.0 with openmpi/4.1.2-gcc-11.2.0 compilers"
  exit 1
elif  [[ "${CLEO_COMPILERNAME}" == "intel" &&
       "${CLEO_CXX_COMPILER}" != "/sw/spack-levante/openmpi-4.1.6-ux3zoj/bin/mpic++" ]]
then
  echo "YAC currently requires intel-oneapi-compilers/2024.2.1-gcc-13.3.0 with\
  openmpi/4.1.6-oneapi-2024.2.1 compilers"
  exit 1
fi
### ---------------------------------------------------- ###

### ------------------ choose YAC build ---------------- ###
source ${bashsrc}/levante_packages.sh
if  [ "${CLEO_COMPILERNAME}" == "gcc" ]
then
  module load ${levante_gcc_netcdf_yac}
  spack load ${levante_gcc_openblas_yac}
elif  [ "${CLEO_COMPILERNAME}" == "intel" ]
then
  module load ${levante_intel_netcdf_yac}
  spack load ${levante_intel_openblas_yac}
fi
export CLEO_YAC_FLAGS="-DCLEO_YAC_MODULE_PATH="${cleo_yac_module_path}" \
  -DCLEO_YAXT_ROOT=${CLEO_YACYAXTROOT}/yaxt \
  -DCLEO_YAC_ROOT=${CLEO_YACYAXTROOT}/yac"
### ---------------------------------------------------- ###
