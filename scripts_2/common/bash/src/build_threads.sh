#!/bin/bash

set -e

COMMON_BASH_SRC="${CLEO_PATH2CLEO}/scripts_2/common/bash/src"
 ### -------------------- check inputs ------------------ ###

 if [ -f "{COMMON_BASH_SRC}/check_inputs.sh" ]; then
   source ${COMMON_BASH_SRC}/check_inputs.sh
 fi
check_args_not_empty "${CLEO_BUILDTYPE}"

if [[ "${CLEO_BUILDTYPE}" != "threads" ]];
then
  echo "Bad inputs, build type for enabling C++ threads on host must be 'threads'"
  exit 1
fi
### --------------------------------------------------- ###

### ------- choose host parallelism kokkos flags ------- ###
export CLEO_KOKKOS_HOST_FLAGS="${CLEO_KOKKOS_HOST_FLAGS} -DKokkos_ENABLE_THREADS=ON"
### ---------------------------------------------------- ###
