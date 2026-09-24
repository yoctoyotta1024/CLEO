#!/bin/bash

set -e
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
bashsrc=${SCRIPT_DIR}

### -------------------- check inputs ------------------ ###
source ${bashsrc}/check_inputs.sh
check_args_not_empty "${CLEO_BUILDTYPE}" "${CLEO_COMPILERNAME} ${CLEO_CXX_COMPILER}"

if [ "${CLEO_BUILDTYPE}" != "cuda" ]
then
  echo "Bad inputs, build type for enabling cuda device must be 'cuda'"
  exit 1
fi

if  [ "${CLEO_COMPILERNAME}" == "intel" ]
then
  echo "CUDA build currently not compatible with intel compiler" # TODO(CB): fix this incompatibility
  exit 1
fi
### ---------------------------------------------------- ###

### --------------- choose CUDA compiler --------------- ###
source ${bashsrc}/levante_packages.sh
spack load ${levante_gcc_cuda}

# set nvcc compiler used by Kokkos nvcc wrapper as CUDA_ROOT/bin/nvcc
cuda_root=${levante_gcc_cuda_root}

# set default (C++) compiler used by kokkos nvcc wrapper
# (wrapper is found in bin directory of Kokkos after its
# installation e.g. build/_deps/kokkos-src/bin/nvcc wrapper)
nvcc_wrapper_default_compiler=${CLEO_CXX_COMPILER}
### ---------------------------------------------------- ###

### ------- choose device parallelism kokkos flags ------- ###
export CLEO_KOKKOS_DEVICE_FLAGS="${CLEO_KOKKOS_DEVICE_FLAGS} -DKokkos_ENABLE_CUDA=ON \
  -DKokkos_ENABLE_CUDA_CONSTEXPR=ON -DKokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE=ON \
  -DCUDA_ROOT=${cuda_root} -DNVCC_WRAPPER_DEFAULT_COMPILER=${nvcc_wrapper_default_compiler}"
### ---------------------------------------------------- ###
