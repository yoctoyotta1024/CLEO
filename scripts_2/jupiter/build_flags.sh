#!/bin/bash

set -e

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)

COMMON_BASH_SRC="${SCRIPT_DIR}/../common"
jupiter_HELPERS_DIR="${SCRIPT_DIR}/helpers"

configure_machine_build_flags() {

    source "${COMMON_BASH_SRC}/check_inputs.sh"
    source "${jupiter_HELPERS_DIR}/jupiter_packages.sh"

    check_args_not_empty \
        "${CLEO_COMPILERNAME}" \
        "${CLEO_ENABLEDEBUG}" \
        "${CLEO_BUILDTYPE}"

    if [[ "${CLEO_COMPILERNAME}" != "gcc" ]]; then
        echo "Error: unsupported compiler '${CLEO_COMPILERNAME}'. Must be 'gcc'."
        exit 1
    fi

    case "${CLEO_ENABLEDEBUG}" in
        true|false) ;;
        *)
            echo "Error: CLEO_ENABLEDEBUG must be 'true' or 'false'."
            exit 1
            ;;
    esac

    jupiter_reset_modules
    jupiter_load_build_stack \
        "${CLEO_COMPILERNAME}" \
        "${CLEO_BUILDTYPE}"

    if ! command -v mpic++ &>/dev/null; then
        echo "Error: 'mpic++' command not found after loading Jupiter toolchain."
        exit 1
    fi

    if ! command -v mpicc &>/dev/null; then
        echo "Error: 'mpicc' command not found after loading Jupiter toolchain."
        exit 1
    fi

    export CLEO_CXX_COMPILER="$(command -v mpic++)"
    export CLEO_CC_COMPILER="$(command -v mpicc)"

    if [[ "${CLEO_ENABLEDEBUG}" == "true" ]]; then

        export CLEO_CXX_FLAGS="${CLEO_CXX_FLAGS} \
            -Werror -Wall -Wextra -pedantic \
            -Wno-unused-parameter -g -gdwarf-4 -O0"

    else

        export CLEO_CXX_FLAGS="${CLEO_CXX_FLAGS} \
            -Werror -Wall -Wextra -pedantic \
            -Wno-unused-parameter -O3"

    fi

    export CLEO_KOKKOS_BASIC_FLAGS="${CLEO_KOKKOS_BASIC_FLAGS} \
        -DKokkos_ARCH_ARMV9_GRACE=ON \
        -DKokkos_ENABLE_SERIAL=ON"

    case "${CLEO_BUILDTYPE}" in

        openmp)
            source "${COMMON_BASH_SRC}/build_openmp.sh"
            configure_openmp_build
            ;;

        threads)
            source "${COMMON_BASH_SRC}/build_threads.sh"
            configure_threads_build
            ;;

        cuda)
            source "${COMMON_BASH_SRC}/build_cuda.sh"
            configure_cuda_build
            export CLEO_KOKKOS_HOST_FLAGS="${CLEO_KOKKOS_HOST_FLAGS} \
                -DKokkos_ENABLE_OPENMP=ON"
            export CLEO_KOKKOS_DEVICE_FLAGS="${CLEO_KOKKOS_DEVICE_FLAGS} \
                -DKokkos_ARCH_HOPPER90=ON"
            ;;

        *)
            echo "Error: unsupported build type '${CLEO_BUILDTYPE}'."
            exit 1
            ;;

    esac
}

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    configure_machine_build_flags "$@"
fi
