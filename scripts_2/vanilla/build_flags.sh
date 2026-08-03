#!/bin/bash

set -e

COMMON_BASH_SRC="${CLEO_PATH2CLEO}/scripts_2/common"

configure_machine_build_flags() {
  ### -------------------- check inputs ------------------ ###
  source ${COMMON_BASH_SRC}/check_inputs.sh
  check_args_not_empty "${CLEO_ENABLEDEBUG}"
  ### ---------------------------------------------------- ###

  ### -------- choose compiler(s) and their flags -------- ###

  if [ "${CLEO_COMPILERNAME}" == "gcc" ]
  then
    # Find a working mpic++/mpicc — skip broken wrappers (e.g. from Anaconda)
    _find_working_mpi() {
      local cmd="$1"
      local found
      # Walk PATH entries in order; test each candidate with --version
      IFS=: read -ra _path_dirs <<< "${PATH}"
      for _dir in "${_path_dirs[@]}"; do
        found="${_dir}/${cmd}"
        if [[ -x "${found}" ]] && "${found}" --version &>/dev/null; then
          echo "${found}"
          return 0
        fi
      done
      return 1
    }

    _mpicxx=$(_find_working_mpi mpic++)
    if [[ -z "${_mpicxx}" ]]; then
      echo "Error: no working 'mpic++' found in PATH. Please ensure a working MPI is installed."
      exit 1
    fi
    _mpicc=$(_find_working_mpi mpicc)
    if [[ -z "${_mpicc}" ]]; then
      echo "Error: no working 'mpicc' found in PATH. Please ensure a working MPI is installed."
      exit 1
    fi

    export CLEO_CXX_COMPILER="${_mpicxx}"
    export CLEO_CC_COMPILER="${_mpicc}"

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
    if declare -F configure_openmp_build >/dev/null; then
      configure_openmp_build
    fi
  elif [[ "${CLEO_BUILDTYPE}" == "threads" ]]; then
    source ${COMMON_BASH_SRC}/build_threads.sh
    if declare -F configure_threads_build >/dev/null; then
      configure_threads_build
    fi
  fi
}

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  configure_machine_build_flags "$@"
fi
