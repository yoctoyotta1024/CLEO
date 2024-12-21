#!/bin/bash

set -e
bashsrc=${CLEO_PATH2CLEO}/scripts/bash/src

### -------------------- check inputs ------------------ ###
source ${bashsrc}/check_inputs.sh
check_args_not_empty "${CLEO_ENABLEYAC}"
### ---------------------------------------------------- ###

### --------------- YAC runtime settings --------------- ###
if [ "${CLEO_ENABLEYAC}" == "true" ]
then
  echo "TODO(CB): something to do with YAC b4 runnnning"
  exit 1
fi
### ---------------------------------------------------- ###
