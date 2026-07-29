#!/bin/bash

set -e

source /etc/profile

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &>/dev/null && pwd )

source "${SCRIPT_DIR}/levante_packages.sh"
source "${CLEO_PATH2CLEO}/scripts_2/common/bash/helpers/build_yac.sh"

case "${CLEO_COMPILERNAME}" in

    gcc)

        if [[ "${CLEO_CXX_COMPILER}" != "/sw/spack-levante/openmpi-4.1.2-mnmady/bin/mpic++" ]]; then
            echo "YAC currently requires gcc/11.2.0 + OpenMPI 4.1.2."
            exit 1
        fi

        module load "${levante_gcc_netcdf_yac}"
        spack load "${levante_gcc_openblas_yac}"

        fyamllib="${levante_gcc_fyamllib}"
        ;;

    intel)

        if [[ "${CLEO_CXX_COMPILER}" != "/sw/spack-levante/openmpi-4.1.6-ux3zoj/bin/mpic++" ]]; then
            echo "YAC currently requires Intel 2024.2.1 + OpenMPI 4.1.6."
            exit 1
        fi

        module load "${levante_intel_netcdf_yac}"
        spack load "${levante_intel_openblas_yac}"

        fyamllib="${levante_intel_fyamllib}"
        ;;

    *)

        echo "Unsupported compiler '${CLEO_COMPILERNAME}'."
        exit 1
        ;;

esac

build_yac "${fyamllib}"
