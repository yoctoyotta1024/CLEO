#!/bin/bash

set -e
bashsrc=${CLEO_PATH2CLEO}/scripts/bash/src

### -------------------- check inputs ------------------ ###
source ${bashsrc}/check_inputs.sh
check_args_not_empty "${CLEO_BUILDTYPE}"
### ---------------------------------------------------- ###

### --------------- set runtime optimisations----------- ###
export OMP_PROC_BIND=spread
export OMP_PLACES=threads

echo "TODO(CB): add runtime optimisations"
exit 1
### ---------------------------------------------------- ###
