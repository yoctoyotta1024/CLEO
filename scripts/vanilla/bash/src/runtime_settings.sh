#!/bin/bash

set -e
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
bashsrc=${SCRIPT_DIR}

### -------------------- check inputs ------------------ ###
source ${bashsrc}/check_inputs.sh
check_args_not_empty "${CLEO_BUILDTYPE}"
check_args_not_empty "${CLEO_COMPILERNAME}" "${CLEO_YACYAXTROOT}"
### ---------------------------------------------------- ###

### --------------- YAC runtime settings --------------- ###
source ${bashsrc}/vanilla_packages.sh
if [ "${CLEO_COMPILERNAME}" == "gcc" ]
then
  fyamllib=${vanilla_gcc_fyamllib}
else
  echo "Bad inputs, YAC on 'vanilla' computer only compatible with "gcc" compiler name"
  exit 1
fi
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${fyamllib}
export PYTHONPATH=${PYTHONPATH}:${CLEO_YACYAXTROOT}/yac/python # path to YAC python bindings
### ---------------------------------------------------- ###


### --------------- set runtime optimisations----------- ###
export OMP_PROC_BIND=spread
export OMP_PLACES=threads

export OMPI_MCA_btl="tcp,self"
export OMPI_MCA_io=ompio
### ---------------------------------------------------- ###
