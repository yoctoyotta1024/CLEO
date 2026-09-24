#!/bin/bash

### Sets the runtime environment for running CLEO on vanilla.
### Usage: configure_machine_runtime_settings [stacksize_limit (kB)]
### Requires CLEO_BUILDTYPE, CLEO_COMPILERNAME and CLEO_YACYAXTROOT to be exported.

set -e

configure_machine_runtime_settings() {
  local stacksize_limit="${1:-}"  # kB

  local machine_dir
  machine_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
  local common_dir="${machine_dir}/../common"

  ### -------------------- check inputs ------------------ ###
  source "${common_dir}/check_inputs.sh"
  check_args_not_empty "${CLEO_BUILDTYPE}" "${CLEO_COMPILERNAME}" "${CLEO_YACYAXTROOT}"
  ### ---------------------------------------------------- ###

  ### ------------ load compiler/runtime stack ------------ ###
  source "${machine_dir}/helpers/vanilla_packages.sh"
  vanilla_load_runtime_stack "${CLEO_COMPILERNAME}" "${CLEO_BUILDTYPE}"
  ### ---------------------------------------------------- ###

  ### --------------- YAC runtime settings --------------- ###
  local fyamllib
  fyamllib=$(vanilla_fyamllib_for_compiler "${CLEO_COMPILERNAME}")
  export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${fyamllib}
  export PYTHONPATH=${PYTHONPATH}:${CLEO_YACYAXTROOT}/yac/python
  ### ---------------------------------------------------- ###

  ### ------------ communication runtime (MPI) ------------ ###
  export OMPI_MCA_btl="tcp,self"
  export OMPI_MCA_io=ompio
  ### ---------------------------------------------------- ###

  ### ------------------ threading ----------------------- ###
  export OMP_PROC_BIND=spread
  export OMP_PLACES=threads
  ### ---------------------------------------------------- ###

  ### ------------------ process limits ------------------ ###
  if [[ -n "${stacksize_limit}" ]]; then
    ulimit -s "${stacksize_limit}"
  fi

  ### ---------------------------------------------------- ###
}

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  configure_machine_runtime_settings "$@"
fi
