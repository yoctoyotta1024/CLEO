#!/bin/bash

set -e

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &>/dev/null && pwd )

configure_machine_yac_flags() {
    source "${SCRIPT_DIR}/vanilla_packages.sh"
    source "${CLEO_PATH2CLEO}/scripts_2/common/build_yac.sh"

    case "${CLEO_COMPILERNAME}" in
        gcc)
            fyamllib="${vanilla_gcc_fyamllib}"
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
