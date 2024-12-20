#!/bin/bash

### -------------------- check inputs ------------------ ###
if [ "${CLEO_BUILDTYPE}" != "" ]
then
  echo "Bad inputs, build type required for runtime optimisations"
  exit 1
fi
### ---------------------------------------------------- ###

### --------------- set runtime optimisations----------- ###
export OMP_PROC_BIND=spread
export OMP_PLACES=threads

echo "TODO(CB): add runtime optimisations"
exit 1
### ---------------------------------------------------- ###
