#!/bin/bash

set -e

### -------------- GCC compiler(s) Packages ------------ ###
jupiter_gcc="GCC/14.3.0" # module load
jupiter_gcc_openmpi="OpenMPI/5.0.8" # module load (loaded after GCC)
jupiter_gcc_cmake="CMake/4.0.3" # module load
jupiter_gcc_cuda="CUDA/13" # module load (auto-loaded as a dependency of GCC on GPU nodes)
jupiter_gcc_cuda_root="/e/software/default/stages/2026/software/CUDA/13" # [cuda_root]/bin/nvcc (match `echo $EBROOTCUDA` after loading CUDA)

### specific gcc compiler compatible packages for YAC installation and usage
jupiter_gcc_netcdf_yac="netCDF/4.9.3" # module load
jupiter_gcc_openblas_yac="OpenBLAS/0.3.30" # module load
jupiter_gcc_fyaml_root="/e/software/default/stages/2026/software/libfyaml/0.9-GCCcore-14.3.0" # no module available, use install path directly (match GCCcore version)
jupiter_gcc_fyamllib="${jupiter_gcc_fyaml_root}/lib"
### specific packages for YAC installation only
jupiter_gcc_netcdf_root="/e/software/default/stages/2026/software/netCDF/4.9.3-gompi-2025b" # match `echo $EBROOTNETCDF` after loading netCDF
### ---------------------------------------------------- ###

### needed so venv pythons linked against libpython can resolve it (module sets LD_LIBRARY_PATH)
jupiter_gcc_python="Python/3.13.5" # module load

jupiter_reset_modules() {
  source /etc/profile
  module --force purge
  module load Stages/2026
  module load StdEnv
}

jupiter_fyamllib_for_compiler() {
  local compilername="$1"
  case "${compilername}" in
    gcc)
      echo "${jupiter_gcc_fyamllib}"
      ;;
    *)
      echo "Error: unsupported compiler '${compilername}'. Must be 'gcc'."
      return 1
      ;;
  esac
}

jupiter_load_build_stack() {
  local compilername="$1"
  local buildtype="$2"

  case "${compilername}" in
    gcc)
      module load "${jupiter_gcc}" "${jupiter_gcc_openmpi}"
      module load "${jupiter_gcc_cmake}"
      if [[ "${buildtype}" == "cuda" ]]; then
        module load "${jupiter_gcc_cuda}"
        export CLEO_CUDA_ROOT="${CLEO_CUDA_ROOT:-${jupiter_gcc_cuda_root}}"
      fi
      ;;
    *)
      echo "Error: unsupported compiler '${compilername}'. Must be 'gcc'."
      return 1
      ;;
  esac
}

jupiter_load_runtime_stack() {
  local compilername="$1"
  local buildtype="$2"

  case "${compilername}" in
    gcc)
      module load "${jupiter_gcc}" "${jupiter_gcc_openmpi}"
      if [[ "${buildtype}" == "cuda" ]]; then
        module load "${jupiter_gcc_cuda}"
        export CLEO_CUDA_ROOT="${CLEO_CUDA_ROOT:-${jupiter_gcc_cuda_root}}"
      fi
      ;;
    *)
      echo "Error: unsupported compiler '${compilername}'. Must be 'gcc'."
      return 1
      ;;
  esac
}

jupiter_load_yac_dependencies() {
  local compilername="$1"

  case "${compilername}" in
    gcc)
      module load "${jupiter_gcc_netcdf_yac}"
      module load "${jupiter_gcc_openblas_yac}"
      ;;
    *)
      echo "Error: unsupported compiler '${compilername}'. Must be 'gcc'."
      return 1
      ;;
  esac
}

jupiter_load_python() {
  local compilername="$1"

  case "${compilername}" in
    gcc)
      module load "${jupiter_gcc_python}"
      ;;
    *)
      echo "Error: unsupported compiler '${compilername}'. Must be 'gcc'."
      return 1
      ;;
  esac
}
