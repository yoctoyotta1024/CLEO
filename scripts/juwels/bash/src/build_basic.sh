#!/bin/bash

set -e
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
bashsrc=${SCRIPT_DIR}

### -------------------- check inputs ------------------ ###
source ${bashsrc}/check_inputs.sh
check_args_not_empty "${CLEO_COMPILERNAME}" "${CLEO_ENABLEDEBUG}"
check_buildtype
check_compilername
### ---------------------------------------------------- ###

### -------- choose compiler(s) and their flags -------- ###
source ${bashsrc}/juwels_packages.sh

if [ "${CLEO_COMPILERNAME}" == "intel" ]
then
  module load ${juwels_intel}
  module load ${juwels_intel_mpi}
  module load ${juwels_intel_cmake}
  export CLEO_CXX_COMPILER=${juwels_icpc_compiler}
  export CLEO_CC_COMPILER=${juwels_icc_compiler}

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
  module load ${juwels_gcc}
  module load ${juwels_gcc_mpi}
  module load ${juwels_gcc_cmake}
  export CLEO_CXX_COMPILER=${juwels_gxx_compiler}
  export CLEO_CC_COMPILER=${juwels_gcc_compiler}

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
