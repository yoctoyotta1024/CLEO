#!/bin/bash

set -e

### -------------- GCC compiler(s) Packages ------------ ###
### specific gcc compiler compatible packages for YAC installation and usage
### hints: https://dkrz-sw.gitlab-pages.dkrz.de/yac/d1/d9f/installing_yac.html
vanilla_gcc_fyaml_root="/opt/homebrew" # match libfyaml
vanilla_gcc_fyamllib="${vanilla_gcc_fyaml_root}/lib"
### specific packages for YAC installation only
vanilla_gcc_netcdf_root="/opt/homebrew" # match netcdf
### ---------------------------------------------------- ###

### A vanilla machine has no module system: the toolchain (mpic++, mpicc,
### cmake) must already be on PATH. These functions keep the same interface
### as the other machines' packages files.

vanilla_reset_modules() {
  :
}

vanilla_fyamllib_for_compiler() {
  local compilername="$1"
  case "${compilername}" in
    gcc)
      echo "${vanilla_gcc_fyamllib}"
      ;;
    *)
      echo "Error: unsupported compiler '${compilername}'. Must be 'gcc'."
      return 1
      ;;
  esac
}

vanilla_load_build_stack() {
  local compilername="$1"
  local buildtype="$2"

  case "${compilername}" in
    gcc)
      if [[ "${buildtype}" == "cuda" ]]; then
        echo "Error: CUDA builds are not supported on a vanilla machine."
        return 1
      fi
      ;;
    *)
      echo "Error: unsupported compiler '${compilername}'. Must be 'gcc'."
      return 1
      ;;
  esac
}

vanilla_load_runtime_stack() {
  vanilla_load_build_stack "$@"
}

vanilla_load_yac_dependencies() {
  local compilername="$1"
  vanilla_fyamllib_for_compiler "${compilername}" > /dev/null
}
