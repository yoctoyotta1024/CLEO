#!/bin/bash

### Sets the YAC cmake flags for building CLEO on jupiter.
### Requires CLEO_COMPILERNAME, CLEO_CXX_COMPILER and CLEO_YACYAXTROOT to be exported.

set -e

configure_machine_yac_flags() {
  local helpers_dir
  helpers_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
  local common_dir="${helpers_dir}/../../common"

  source "${helpers_dir}/jupiter_packages.sh"
  source "${common_dir}/build_yac.sh"
  jupiter_load_yac_dependencies "${CLEO_COMPILERNAME}"

  ### ---- check compiler is compatible with YAC install ---- ###
  case "${CLEO_COMPILERNAME}" in
    gcc)
      if [[ "${CLEO_CXX_COMPILER}" != "/e/software/default/stages/2026/software/OpenMPI/5.0.8-GCC-14.3.0/bin/mpic++" ]]; then
        echo "YAC currently requires GCC/14.3.0 + OpenMPI 5.0.8."
        exit 1
      fi
      ;;
    *)
      echo "Unsupported compiler '${CLEO_COMPILERNAME}' for YAC on jupiter."
      exit 1
      ;;
  esac
  ### ------------------------------------------------------ ###

  local fyamllib
  fyamllib=$(jupiter_fyamllib_for_compiler "${CLEO_COMPILERNAME}")
  build_yac "${fyamllib}"
}

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  configure_machine_yac_flags "$@"
fi
