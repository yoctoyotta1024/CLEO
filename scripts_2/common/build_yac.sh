#!/bin/bash

set -e

build_yac() {

    local fyamllib="$1"

    local common_dir="${CLEO_PATH2CLEO}/scripts_2/common"
    local cleo_yac_module_path="${CLEO_PATH2CLEO}/libs/coupldyn_yac/cmake"

    source "${common_dir}/check_inputs.sh"

    check_args_not_empty \
        "${CLEO_YACYAXTROOT}" \
        "${fyamllib}"

    export CLEO_YAC_FLAGS="
        -DCLEO_YAC_MODULE_PATH=${cleo_yac_module_path}
        -DCLEO_FYAMLLIB=${fyamllib}
        -DCLEO_YAXT_ROOT=${CLEO_YACYAXTROOT}/yaxt
        -DCLEO_YAC_ROOT=${CLEO_YACYAXTROOT}/yac"
}
