#!/bin/bash

set -e
source /etc/profile
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
bashsrc=${SCRIPT_DIR}

### -------------------- check inputs ------------------ ###
source ${bashsrc}/check_inputs.sh
check_args_not_empty "${CLEO_COMPILERNAME}" "${CLEO_ENABLEDEBUG}"
check_buildtype
check_compilername
### ---------------------------------------------------- ###

### -------- choose compiler(s) and their flags -------- ###
source ${bashsrc}/levante_packages.sh

if [ "${CLEO_COMPILERNAME}" == "intel" ]
then
  module load ${levante_intel} ${levante_intel_openmpi}
  spack load ${levante_intel_cmake}
  export CLEO_CXX_COMPILER=${levante_icpc_compiler}
  export CLEO_CC_COMPILER=${levante_icc_compiler}

  if [ "${CLEO_ENABLEDEBUG}" == "true" ]
  then
    ### for correctness and debugging (note -gdwarf-4 not possible for nvc++) use:
    export CLEO_CXX_FLAGS="${CLEO_CXX_FLAGS} -Werror -Wall -Wextra \
      -pedantic -Wno-unused-parameter -g -gdwarf-4 -O0" # correctness and debugging
  else
    ### for performance use:
    export CLEO_CXX_FLAGS="${CLEO_CXX_FLAGS} -Werror -Wall -Wextra \
      -pedantic -Wno-unused-parameter -O3 -fma"
  fi
elif [ "${CLEO_COMPILERNAME}" == "gcc" ]
then
  module load ${levante_gcc} ${levante_gcc_openmpi}
  spack load ${levante_gcc_cmake}
  export CLEO_CXX_COMPILER=${levante_gxx_compiler}
  export CLEO_CC_COMPILER=${levante_gcc_compiler}

  if [ "${CLEO_ENABLEDEBUG}" == "true" ]
  then
    ### for correctness and debugging (note -gdwarf-4 not possible for nvc++) use:
    export CLEO_CXX_FLAGS="${CLEO_CXX_FLAGS} -Werror -Wno-unused-parameter -Wall -Wextra \
      -pedantic -g -gdwarf-4 -O0 -mpc64"
  else
    ### for performance use:
    export CLEO_CXX_FLAGS="${CLEO_CXX_FLAGS} -Werror -Wall -Wextra \
      -pedantic -Wno-unused-parameter -O3 -mfma"
  fi
fi
### ---------------------------------------------------- ###

### ------------ choose basic kokkos flags ------------- ###
export CLEO_KOKKOS_BASIC_FLAGS="${CLEO_KOKKOS_BASIC_FLAGS} \
  -DKokkos_ARCH_NATIVE=ON -DKokkos_ENABLE_SERIAL=ON"
### ---------------------------------------------------- ###
