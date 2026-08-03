#!/bin/bash

set -e

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &>/dev/null && pwd )
LEVANTE_DIR="${SCRIPT_DIR}/.."
COMMON_DIR="${SCRIPT_DIR}/../../common"

source "${LEVANTE_DIR}/helpers/levante_packages.sh"
levante_reset_modules

path2toolsrepo="$1"
compilername="${2:-gcc}"
root4tools="${3:-${CLEO_KOKKOSTOOLS}}"

case "${compilername}" in

    gcc)
        levante_load_build_stack "${compilername}" "openmp"
        ;;

    intel)
        levante_load_build_stack "${compilername}" "openmp"
        ;;

    *)
        echo "Unknown compiler '${compilername}'."
        exit 1
        ;;
esac

source "${COMMON_DIR}/install_kokkos_tools.sh"

install_kokkos_tools \
    "${path2toolsrepo}" \
    "${root4tools}"
