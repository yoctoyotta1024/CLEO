#!/bin/bash

set -e

COMMON_BASH_SRC="${CLEO_PATH2CLEO}/scripts_2/common/bash/src"
LEVANTE_BASH_SRC="${CLEO_PATH2CLEO}/scripts_2/levante/bash/src"

### -------------------- check inputs ------------------ ###
source ${COMMON_BASH_SRC}/check_inputs.sh
check_args_not_empty "${CLEO_COMPILERNAME}" "${CLEO_ENABLEDEBUG}"
### ---------------------------------------------------- ###

### -------- choose compiler(s) and their flags -------- ###
source ${LEVANTE_BASH_SRC}/levante_packages.sh

if [ "${CLEO_COMPILERNAME}" == "intel" ]
then
  if [ "${CLEO_BUILDTYPE}" == "cuda" ]
  then
    echo "Error: CUDA build on Levante is not compatible with the intel compiler. Use gcc."
    exit 1
  fi
  module load ${levante_intel} ${levante_intel_openmpi}
  spack load ${levante_intel_cmake}
  export CLEO_CXX_COMPILER="$(command -v mpic++)"
  export CLEO_CC_COMPILER="$(command -v mpicc)"

  if [ "${CLEO_ENABLEDEBUG}" == "true" ]
  then
    export CLEO_CXX_FLAGS="${CLEO_CXX_FLAGS} -Werror -Wall -Wextra \
      -pedantic -Wno-unused-parameter -g -gdwarf-4 -O0"
  else
    export CLEO_CXX_FLAGS="${CLEO_CXX_FLAGS} -Werror -Wall -Wextra \
      -pedantic -Wno-unused-parameter -O3 -fma"
  fi

elif [ "${CLEO_COMPILERNAME}" == "gcc" ]
then
  module load ${levante_gcc} ${levante_gcc_openmpi}
  spack load ${levante_gcc_cmake}
  export CLEO_CXX_COMPILER="$(command -v mpic++)"
  export CLEO_CC_COMPILER="$(command -v mpicc)"

  if [ "${CLEO_BUILDTYPE}" == "cuda" ]
  then
    spack load ${levante_gcc_cuda}
    export CLEO_CUDA_ROOT=${levante_gcc_cuda_root}
  fi

  if [ "${CLEO_ENABLEDEBUG}" == "true" ]
  then
    export CLEO_CXX_FLAGS="${CLEO_CXX_FLAGS} -Werror -Wno-unused-parameter -Wall -Wextra \
      -pedantic -g -gdwarf-4 -O0 -mpc64"
  else
    export CLEO_CXX_FLAGS="${CLEO_CXX_FLAGS} -Werror -Wall -Wextra \
      -pedantic -Wno-unused-parameter -O3 -mfma"
  fi

else
  echo "Error: unsupported compiler '${CLEO_COMPILERNAME}'. Must be 'gcc' or 'intel'."
  exit 1
fi
### ---------------------------------------------------- ###

### ------------ choose basic kokkos flags ------------- ###
export CLEO_KOKKOS_BASIC_FLAGS="${CLEO_KOKKOS_BASIC_FLAGS} \
  -DKokkos_ARCH_NATIVE=ON -DKokkos_ENABLE_SERIAL=ON"
### ---------------------------------------------------- ###

### ------ choose host/device parallelism flags --------- ###
if [[ "${CLEO_BUILDTYPE}" == "openmp" ]]; then
  source ${COMMON_BASH_SRC}/build_openmp.sh
elif [[ "${CLEO_BUILDTYPE}" == "threads" ]]; then
  source ${COMMON_BASH_SRC}/build_threads.sh
elif [[ "${CLEO_BUILDTYPE}" == "cuda" ]]; then
  source ${COMMON_BASH_SRC}/build_cuda.sh
fi
### ---------------------------------------------------- ###
