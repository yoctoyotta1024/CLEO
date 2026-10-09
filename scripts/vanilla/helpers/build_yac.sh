#!/bin/bash

### Sets the YAC cmake flags for building CLEO on vanilla.
### Requires CLEO_COMPILERNAME, CLEO_CXX_COMPILER and CLEO_YACYAXTROOT to be exported.

set -e

configure_machine_yac_flags() {
  local helpers_dir
  helpers_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
  local common_dir="${helpers_dir}/../../common"

  source "${helpers_dir}/vanilla_packages.sh"
  source "${common_dir}/build_yac.sh"
  vanilla_load_yac_dependencies "${CLEO_COMPILERNAME}"

  ### ---- check compiler is compatible with YAC install ---- ###
  case "${CLEO_COMPILERNAME}" in
    gcc)
      # no fixed toolchain on a vanilla machine: the YAC install must
      # have been built with the same mpicc/mpic++ found on PATH
      ;;
    *)
      echo "Unsupported compiler '${CLEO_COMPILERNAME}' for YAC on vanilla."
      exit 1
      ;;
  esac
  ### ------------------------------------------------------ ###

  local fyamllib
  fyamllib=$(vanilla_fyamllib_for_compiler "${CLEO_COMPILERNAME}")
  build_yac "${fyamllib}"
}

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  configure_machine_yac_flags "$@"
fi
