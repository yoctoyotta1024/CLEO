#!/bin/bash

set -e

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &>/dev/null && pwd )
COMMON_DIR="${SCRIPT_DIR}/../../common"

configure_machine_yac_flags() {
    source /etc/profile

    source "${SCRIPT_DIR}/levante_packages.sh"
    source "${COMMON_DIR}/build_yac.sh"
    levante_load_yac_dependencies "${CLEO_COMPILERNAME}"

    case "${CLEO_COMPILERNAME}" in

        gcc)

            if [[ "${CLEO_CXX_COMPILER}" != "/sw/spack-levante/openmpi-4.1.2-mnmady/bin/mpic++" ]]; then
                echo "YAC currently requires gcc/11.2.0 + OpenMPI 4.1.2."
                exit 1
            fi

            fyamllib=$(levante_fyamllib_for_compiler "${CLEO_COMPILERNAME}")
            ;;

        intel)

            if [[ "${CLEO_CXX_COMPILER}" != "/sw/spack-levante/openmpi-4.1.6-ux3zoj/bin/mpic++" ]]; then
                echo "YAC currently requires Intel 2024.2.1 + OpenMPI 4.1.6."
                exit 1
            fi

            fyamllib=$(levante_fyamllib_for_compiler "${CLEO_COMPILERNAME}")
            ;;

        *)

            echo "Unsupported compiler '${CLEO_COMPILERNAME}'."
            exit 1
            ;;

    esac

    build_yac "${fyamllib}"
}

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    configure_machine_yac_flags "$@"
fi
