#!/bin/bash

set -e

### -------------------- check inputs ------------------ ###
if [ "${CLEO_BUILDTYPE}" != "serial" ] &&
   [ "${CLEO_BUILDTYPE}" != "openmp" ] &&
   [ "${CLEO_BUILDTYPE}" != "threads" ] &&
   [ "${CLEO_BUILDTYPE}" != "cuda" ];
then
  echo "Bad inputs, build type must be 'serial', 'openmp', 'threads' or 'cuda'"
  exit 1
fi

if [ "${CLEO_COMPILERNAME}" == "" || "${CLEO_ENABLEDEBUG}" ]
then
  echo "Bad inputs, please specify CLEO compiler name and debuging mode"
  exit 1
fi

if [ "${CLEO_COMPILERNAME}" != "intel" ] &&
   [ "${CLEO_COMPILERNAME}" != "gcc" ];
then
  echo "Bad inputs, CLEO compiler name must be 'intel' or 'gcc'"
  exit 1
fi
### ---------------------------------------------------- ###

### -------- choose compiler(s) and their flags -------- ###
if [ "${CLEO_COMPILERNAME}" == "intel" ]
then
  echo "TODO(CB): intel compiler support"
  exit 1
  if [ "${CLEO_ENABLEDEBUG}" == "true" ]
  then
    ### for correctness and debugging (note -gdwarf-4 not possible for nvc++) use:
    export CLEO_CXX_FLAGS="${CLEO_CXX_FLAGS} -Werror -Wno-unused-parameter -Wall -Wextra -pedantic -g -gdwarf-4 -O0"
  else
    ### for performance use:
    export CLEO_CXX_FLAGS="${CLEO_CXX_FLAGS} -Werror -Wall -pedantic -O3"
  fi
fi

elif [ "${CLEO_COMPILERNAME}" == "gcc" ]
then
  module load gcc/11.2.0-gcc-11.2.0 openmpi/4.1.2-gcc-11.2.0
  spack load cmake@3.23.1%gcc
  export CLEO_CXX_COMPILER="/sw/spack-levante/openmpi-4.1.2-mnmady/bin/mpic++"
  export CLEO_CC_COMPILER="/sw/spack-levante/openmpi-4.1.2-mnmady/bin/mpicc"

  echo "TODO(CB): update gcc compiler version (in YAC and cuda too!)"
  exit 1

  if [ "${CLEO_ENABLEDEBUG}" == "true" ]
  then
    ### for correctness and debugging (note -gdwarf-4 not possible for nvc++) use:
    export CLEO_CXX_FLAGS="${CLEO_CXX_FLAGS} -Werror -Wno-unused-parameter -Wall -Wextra -pedantic -g -gdwarf-4 -O0 -mpc64"
  else
    ### for performance use:
    export CLEO_CXX_FLAGS="${CLEO_CXX_FLAGS} -Werror -Wall -pedantic -O3"
  fi
fi
### ---------------------------------------------------- ###

### ------------ choose basic kokkos flags ------------- ###
export CLEO_KOKKOS_BASIC_FLAGS="${CLEO_KOKKOS_BASIC_FLAGS} -DKokkos_ARCH_NATIVE=ON -DKokkos_ENABLE_SERIAL=ON"
### ---------------------------------------------------- ###
