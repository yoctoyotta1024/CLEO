#!/bin/bash

### -------------------- check inputs ------------------ ###
if [ "${CLEO_BUILDTYPE}" != "threads" ] &&
   [ "${CLEO_BUILDTYPE}" != "cuda" ];
then
  echo "Bad inputs, build type for enabling C++ threads on host must be 'threads' or 'cuda'"
  exit 1
fi
### --------------------------------------------------- ###

### ------- choose host parallelism kokkos flags ------- ###
export CLEO_KOKKOS_HOST_FLAGS="${CLEO_KOKKOS_HOST_FLAGS} -DKokkos_ENABLE_THREADS=ON"
### ---------------------------------------------------- ###
