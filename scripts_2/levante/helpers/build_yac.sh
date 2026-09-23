#!/bin/bash

### Sets the YAC cmake flags for building CLEO on levante.
### Requires CLEO_COMPILERNAME, CLEO_CXX_COMPILER and CLEO_YACYAXTROOT to be exported.

set -e

configure_machine_yac_flags() {
  local helpers_dir
  helpers_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
  local common_dir="${helpers_dir}/../../common"

  source "${helpers_dir}/levante_packages.sh"
  source "${common_dir}/build_yac.sh"
  levante_load_yac_dependencies "${CLEO_COMPILERNAME}"

  ### ---- check compiler is compatible with YAC install ---- ###
  case "${CLEO_COMPILERNAME}" in
    gcc)
      if [[ "${CLEO_CXX_COMPILER}" != "/sw/spack-levante/openmpi-4.1.2-mnmady/bin/mpic++" ]]; then
        echo "YAC currently requires gcc/11.2.0 + OpenMPI 4.1.2."
        exit 1
      fi
      ;;
    intel)
      if [[ "${CLEO_CXX_COMPILER}" != "/sw/spack-levante/openmpi-4.1.6-ux3zoj/bin/mpic++" ]]; then
        echo "YAC currently requires Intel 2024.2.1 + OpenMPI 4.1.6."
        exit 1
      fi
      ;;
    *)
      echo "Unsupported compiler '${CLEO_COMPILERNAME}' for YAC on levante."
      exit 1
      ;;
  esac
  ### ------------------------------------------------------ ###

  local fyamllib
  fyamllib=$(levante_fyamllib_for_compiler "${CLEO_COMPILERNAME}")
  build_yac "${fyamllib}"
}

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  configure_machine_yac_flags "$@"
fi
