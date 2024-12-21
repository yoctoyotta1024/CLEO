#!/bin/bash

set -e
bashsrc=${CLEO_PATH2CLEO}/scripts/bash/src

### -------------------- check inputs ------------------ ###
source ${bashsrc}/check_inputs.sh
check_args_not_empty "${CLEO_BUILDTYPE}" "${CLEO_CXX_COMPILER}"

if [ "${CLEO_BUILDTYPE}" != "cuda" ]
then
  echo "Bad inputs, build type for enabling cuda device must be 'cuda'"
  exit 1
fi

if  [ "${CLEO_CXX_COMPILER}" == "intel" ]
then
  echo "CUDA build currently not compatible with intel compiler" # TODO(CB): fix this incompatibility
  exit 1
fi
### ---------------------------------------------------- ###

### --------------- choose CUDA compiler --------------- ###
# set nvcc compiler used by Kokkos nvcc wrapper as CUDA_ROOT/bin/nvcc
# NOTE(!) this path should correspond to the loaded nvhpc module.
# you can get a clue for the correct path e.g. via 'spack find -p nvhpc@23.9'
cuda_root="/sw/spack-levante/nvhpc-23.9-xpxqeo/Linux_x86_64/23.9/cuda/"

# set default (C++) compiler used by kokkos nvcc wrapper
# (wrapper is found in bin directory of Kokkos after its
# installation e.g. build/_deps/kokkos-src/bin/nvcc wrapper)
nvcc_wrapper_default_compiler=${CLEO_CXX_COMPILER}
### ---------------------------------------------------- ###

### ------- choose device parallelism kokkos flags ------- ###
module load nvhpc/23.9-gcc-11.2.0 # must match gcc compiler (!)
export CLEO_KOKKOS_DEVICE_FLAGS="${CLEO_KOKKOS_DEVICE_FLAGS} -DKokkos_ENABLE_CUDA=ON \
  -DKokkos_ENABLE_CUDA_CONSTEXPR=ON -DKokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE=ON \
  -DCUDA_ROOT=${cuda_root} -DNVCC_WRAPPER_DEFAULT_COMPILER=${nvcc_wrapper_default_compiler}"
### ---------------------------------------------------- ###
