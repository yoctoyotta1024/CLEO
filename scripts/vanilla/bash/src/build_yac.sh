#!/bin/bash

set -e
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
bashsrc=${SCRIPT_DIR}
cleo_yac_module_path="${CLEO_PATH2CLEO}/libs/coupldyn_yac/cmake"

### -------------------- check inputs ------------------ ###
source ${bashsrc}/check_inputs.sh
check_args_not_empty "${CLEO_PATH2CLEO}" "${CLEO_COMPILERNAME}" "${CLEO_CXX_COMPILER}" \
  "${CLEO_YACYAXTROOT}" "${CLEO_FYAMLLIB}"
check_compilername
### ---------------------------------------------------- ###

### ------------------ choose YAC build ---------------- ###
export CLEO_YAC_FLAGS="-DCLEO_YAC_MODULE_PATH="${cleo_yac_module_path}" \
  -DCLEO_FYAMLLIB=${CLEO_FYAMLLIB} \
  -DCLEO_YAXT_ROOT=${CLEO_YACYAXTROOT}/yaxt \
  -DCLEO_YAC_ROOT=${CLEO_YACYAXTROOT}/yac"
### ---------------------------------------------------- ###
