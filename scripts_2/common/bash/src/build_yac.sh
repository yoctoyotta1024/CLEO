#!/bin/bash

set -e

COMMON_BASH_SRC="${CLEO_PATH2CLEO}/scripts_2/common/bash/src"
cleo_yac_module_path="${CLEO_PATH2CLEO}/libs/coupldyn_yac/cmake"

### -------------------- check inputs ------------------ ###

if [ -f "{COMMON_BASH_SRC}/check_inputs.sh" ]; then
  source ${COMMON_BASH_SRC}/check_inputs.sh
fi
check_args_not_empty "${CLEO_YACYAXTROOT}" "${CLEO_FYAMLLIB}"
### ---------------------------------------------------- ###

### ------------------ choose YAC build ---------------- ###
export CLEO_YAC_FLAGS="-DCLEO_YAC_MODULE_PATH="${cleo_yac_module_path}" \
  -DCLEO_FYAMLLIB=${CLEO_FYAMLLIB} \
  -DCLEO_YAXT_ROOT=${CLEO_YACYAXTROOT}/yaxt \
  -DCLEO_YAC_ROOT=${CLEO_YACYAXTROOT}/yac"
### ---------------------------------------------------- ###
