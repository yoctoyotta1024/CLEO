#!/bin/bash

set -e

COMMON_BASH_SRC="${CLEO_PATH2CLEO}/scripts_2/common/bash/src"

### -------------------- check inputs ------------------ ###
source ${COMMON_BASH_SRC}/check_inputs.sh
check_args_not_empty "${CLEO_BUILDTYPE}" "${CLEO_COMPILERNAME}" \
  "${CLEO_CXX_COMPILER}" "${CLEO_CUDA_ROOT}"

if [ "${CLEO_BUILDTYPE}" != "cuda" ]
then
  echo "Bad inputs, build type for enabling cuda device must be 'cuda'"
  exit 1
fi
### ---------------------------------------------------- ###

### ------- choose CUDA compiler and device flags ------- ###

# set nvcc compiler used by Kokkos nvcc wrapper as CUDA_ROOT/bin/nvcc
cuda_root=${CLEO_CUDA_ROOT}

# set default (C++) compiler used by the kokkos nvcc wrapper
# (wrapper is found in bin directory of Kokkos after its installation,
#  e.g. build/_deps/kokkos-src/bin/nvcc_wrapper)
nvcc_wrapper_default_compiler=${CLEO_CXX_COMPILER}

export CLEO_KOKKOS_DEVICE_FLAGS="${CLEO_KOKKOS_DEVICE_FLAGS} \
  -DKokkos_ENABLE_CUDA=ON \
  -DKokkos_ENABLE_CUDA_CONSTEXPR=ON \
  -DKokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE=ON \
  -DCUDA_ROOT=${cuda_root} \
  -DNVCC_WRAPPER_DEFAULT_COMPILER=${nvcc_wrapper_default_compiler}"
### ---------------------------------------------------- ###
