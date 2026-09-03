#!/bin/bash

set -e

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &>/dev/null && pwd )
COMMON_DIR="${SCRIPT_DIR}/../../common"

configure_machine_yac_flags() {
    source /etc/profile

    source "${SCRIPT_DIR}/jupiter_packages.sh"
    source "${COMMON_DIR}/build_yac.sh"
    jupiter_load_yac_dependencies "${CLEO_COMPILERNAME}"

    case "${CLEO_COMPILERNAME}" in

        gcc)

            if [[ "${CLEO_CXX_COMPILER}" != "/e/software/default/stages/2026/software/OpenMPI/5.0.8-GCC-14.3.0/bin/mpic++" ]]; then
                echo "YAC currently requires GCC/14.3.0 + OpenMPI 5.0.8."
                exit 1
            fi

            fyamllib=$(jupiter_fyamllib_for_compiler "${CLEO_COMPILERNAME}")
            ;;

        *)

            echo "Unsupported compiler '${CLEO_COMPILERNAME}'. Must be 'gcc'."
            exit 1
            ;;

    esac

    build_yac "${fyamllib}"
}

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    configure_machine_yac_flags "$@"
fi
