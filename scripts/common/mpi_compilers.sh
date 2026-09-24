#!/bin/bash

### Selects the MPI compiler wrappers used to build CLEO.
### Must be called AFTER the machine's toolchain (modules) has been loaded.
### Walks PATH in order and takes the first mpic++/mpicc that actually runs,
### skipping broken wrappers (e.g. from Anaconda).

find_working_command() {
  local cmd="$1"
  local dir
  local candidate
  local path_dirs
  IFS=: read -r -a path_dirs <<< "${PATH}"
  for dir in "${path_dirs[@]}"; do
    candidate="${dir}/${cmd}"
    if [[ -x "${candidate}" ]] && "${candidate}" --version &>/dev/null; then
      echo "${candidate}"
      return 0
    fi
  done
  return 1
}

configure_mpi_compilers() {
  local mpicxx
  local mpicc

  mpicxx=$(find_working_command mpic++) || true
  if [[ -z "${mpicxx}" ]]; then
    echo "Error: no working 'mpic++' found in PATH after loading the ${CLEO_MACHINE} toolchain."
    exit 1
  fi

  mpicc=$(find_working_command mpicc) || true
  if [[ -z "${mpicc}" ]]; then
    echo "Error: no working 'mpicc' found in PATH after loading the ${CLEO_MACHINE} toolchain."
    exit 1
  fi

  export CLEO_CXX_COMPILER="${mpicxx}"
  export CLEO_CC_COMPILER="${mpicc}"
}

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  configure_mpi_compilers "$@"
fi
