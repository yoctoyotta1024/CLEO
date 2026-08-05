#!/bin/bash

set -e
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
COMMON_DIR="${SCRIPT_DIR}/../common"
HELPER_DIR="${SCRIPT_DIR}/helpers"

configure_machine_runtime_settings() {
  ### -------------------- check inputs ------------------ ###
  source "${COMMON_DIR}/check_inputs.sh"
  check_args_not_empty "${CLEO_BUILDTYPE}"
  check_args_not_empty "${CLEO_COMPILERNAME}" "${CLEO_YACYAXTROOT}"
  ### ---------------------------------------------------- ###

  ### --------------- YAC runtime settings --------------- ###
  source "${HELPER_DIR}/vanilla_packages.sh"
  if [ "${CLEO_COMPILERNAME}" == "gcc" ]; then
    fyamllib=${vanilla_gcc_fyamllib}
  else
    echo "Bad inputs, YAC on 'vanilla' computer only compatible with 'gcc' compiler name"
    exit 1
  fi
  export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${fyamllib}
  export PYTHONPATH=${PYTHONPATH}:${CLEO_YACYAXTROOT}/yac/python
  ### ---------------------------------------------------- ###


  ### --------------- set runtime optimisations----------- ###
  export OMP_PROC_BIND=spread
  export OMP_PLACES=threads

  export OMPI_MCA_btl="tcp,self"
  export OMPI_MCA_io=ompio
  ### ---------------------------------------------------- ###
}

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  configure_machine_runtime_settings "$@"
fi
