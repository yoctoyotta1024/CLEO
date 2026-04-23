#!/bin/bash

set -e

COMMON_BASH_SRC="${CLEO_PATH2CLEO}/scripts_2/common/bash/src"
### -------------------- check inputs ------------------ ###
source ${COMMON_BASH_SRC}/check_inputs.sh
check_args_not_empty "${CLEO_ENABLEDEBUG}"
### ---------------------------------------------------- ###

### -------- choose compiler(s) and their flags -------- ###

if [ "${CLEO_COMPILERNAME}" == "gcc" ]
then
  if ! command -v mpic++ &> /dev/null; then
    echo "Error: 'mpic++' command not found. Please ensure MPI is installed and available in PATH."
    exit 1
  fi
  if ! command -v mpicc &> /dev/null; then
    echo "Error: 'mpicc' command not found. Please ensure MPI is installed and available in PATH."
    exit 1
  fi

  export CLEO_CXX_COMPILER="$(command -v mpic++)"
  export CLEO_CC_COMPILER="$(command -v mpicc)"

  if [ "${CLEO_ENABLEDEBUG}" == "true" ]
  then
    export CLEO_CXX_FLAGS="${CLEO_CXX_FLAGS} -Werror -Wno-unused-parameter -Wall -Wextra \
      -pedantic -g -gdwarf-4 -O0 -mp64"
  else
    export CLEO_CXX_FLAGS="${CLEO_CXX_FLAGS} -Werror -Wall -Wextra \
      -pedantic -Wno-unused-parameter -O3" # -mfma" # (mfma not compatible with apple silicon arch)
  fi
else
  echo "Error: unsupported compiler '${CLEO_COMPILERNAME}'. Only 'gcc' is supported in this script."
  exit 1
fi

### ---------------------------------------------------- ###

### ------------ choose basic kokkos flags ------------- ###
export CLEO_KOKKOS_BASIC_FLAGS="${CLEO_KOKKOS_BASIC_FLAGS} \
  -DKokkos_ARCH_NATIVE=ON -DKokkos_ENABLE_SERIAL=ON"
### ---------------------------------------------------- ###

if [[ "${CLEO_BUILDTYPE}" == "openmp" ]]; then
  source ${COMMON_BASH_SRC}/build_openmp.sh
elif [[ "${CLEO_BUILDTYPE}" == "threads" ]]; then
  source ${COMMON_BASH_SRC}/build_threads.sh
fi
