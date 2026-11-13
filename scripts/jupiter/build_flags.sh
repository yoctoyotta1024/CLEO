#!/bin/bash

### Sets the compiler, compiler flags and Kokkos flags for building CLEO on jupiter.
### Requires CLEO_COMPILERNAME, CLEO_BUILDTYPE and CLEO_ENABLEDEBUG to be exported.

set -e

configure_machine_build_flags() {
  local machine_dir
  machine_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
  local common_dir="${machine_dir}/../common"

  ### -------------------- check inputs ------------------ ###
  source "${common_dir}/check_inputs.sh"
  check_args_not_empty "${CLEO_COMPILERNAME}" "${CLEO_ENABLEDEBUG}" "${CLEO_BUILDTYPE}"
  check_value_in_list CLEO_ENABLEDEBUG "${CLEO_ENABLEDEBUG}" true false
  ### ---------------------------------------------------- ###

  ### ----- load toolchain and choose MPI compilers ------ ###
  source "${machine_dir}/helpers/jupiter_packages.sh"
  jupiter_reset_modules
  jupiter_load_build_stack "${CLEO_COMPILERNAME}" "${CLEO_BUILDTYPE}"
  jupiter_load_python "${CLEO_COMPILERNAME}" # so a venv python linked against libpython runs (e.g. for CMake)

  source "${common_dir}/mpi_compilers.sh"
  configure_mpi_compilers
  ### ---------------------------------------------------- ###

  ### --------------- choose compiler flags -------------- ###
  case "${CLEO_COMPILERNAME}" in
    gcc)
      if [[ "${CLEO_ENABLEDEBUG}" == "true" ]]; then
        export CLEO_CXX_FLAGS="${CLEO_CXX_FLAGS} -Werror -Wall -Wextra -pedantic -Wno-unused-parameter -g -gdwarf-4 -O0"
      else
        export CLEO_CXX_FLAGS="${CLEO_CXX_FLAGS} -Werror -Wall -Wextra -pedantic -Wno-unused-parameter -O3"
      fi
      ;;
    *)
      echo "Error: unsupported compiler '${CLEO_COMPILERNAME}' on jupiter."
      exit 1
      ;;
  esac
  ### ---------------------------------------------------- ###

  ### ------------ choose basic kokkos flags ------------- ###
  export CLEO_KOKKOS_BASIC_FLAGS="${CLEO_KOKKOS_BASIC_FLAGS} \
    -DKokkos_ARCH_ARMV9_GRACE=ON -DKokkos_ENABLE_SERIAL=ON"
  ### ---------------------------------------------------- ###

  ### ------ choose host/device parallelism flags -------- ###
  case "${CLEO_BUILDTYPE}" in
    serial)
      ;;
    openmp)
      source "${common_dir}/build_openmp.sh"
      configure_openmp_build
      ;;
    threads)
      source "${common_dir}/build_threads.sh"
      configure_threads_build
      ;;
    cuda)
      source "${common_dir}/build_cuda.sh"
      configure_cuda_build
      export CLEO_KOKKOS_HOST_FLAGS="${CLEO_KOKKOS_HOST_FLAGS} \
        -DKokkos_ENABLE_OPENMP=ON"
      export CLEO_KOKKOS_DEVICE_FLAGS="${CLEO_KOKKOS_DEVICE_FLAGS} \
        -DKokkos_ARCH_HOPPER90=ON"
      ;;
    *)
      echo "Error: unsupported build type '${CLEO_BUILDTYPE}' on jupiter."
      exit 1
      ;;
  esac
  ### ---------------------------------------------------- ###
}

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  configure_machine_build_flags "$@"
fi
