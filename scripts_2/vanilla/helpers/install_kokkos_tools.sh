#!/bin/bash

set -e

path2toolsrepo="$1"
compilername="${2:-gcc}"
root4tools="${3:-${CLEO_KOKKOSTOOLS}}"

if [[ "${compilername}" != "gcc" ]]; then
    echo "ERROR: Vanilla currently supports only gcc."
    exit 1
fi

COMMON_DIR="${CLEO_PATH2CLEO}/scripts_2/common"

source "${COMMON_DIR}/install_kokkos_tools.sh"

install_kokkos_tools \
    "${path2toolsrepo}" \
    "${root4tools}"
