#!/bin/bash

set -e

source /etc/profile

module purge
spack unload --all

LEVANTE_DIR="${CLEO_PATH2CLEO}/scripts_2/levante"
COMMON_DIR="${CLEO_PATH2CLEO}/scripts_2/common"

source "${LEVANTE_DIR}/helpers/levante_packages.sh"

path2toolsrepo="$1"
compilername="${2:-gcc}"
root4tools="${3:-${CLEO_KOKKOSTOOLS}}"

case "${compilername}" in

    gcc)
        module load "${levante_gcc}" "${levante_gcc_openmpi}"
        spack load "${levante_gcc_cmake}"
        ;;

    intel)
        module load "${levante_intel}" "${levante_intel_openmpi}"
        spack load "${levante_intel_cmake}"
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
