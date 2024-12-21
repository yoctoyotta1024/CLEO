#!/bin/bash

levante_gcc=gcc/11.2.0-gcc-11.2.0 # bcn7mbu
levante_gcc_cmake=cmake@3.26.3%gcc@=11.2.0/fuvwuhz
levante_gcc_openmpi=openmpi@4.1.2%gcc@11.2.0
levante_gxx_compiler="/sw/spack-levante/openmpi-4.1.2-mnmady/bin/mpic++"
levante_gcc_compiler="/sw/spack-levante/openmpi-4.1.2-mnmady/bin/mpicc"

levante_gcc_cuda=cuda@12.2.0%gcc@=11.2.0
levante_gcc_cuda_root="/sw/spack-levante/cuda-12.2.0-2ttufp/"
# NOTE(!) this path should correspond to the loaded nvhpc module.
# you can get a clue for the correct path e.g. via 'spack find -p nvhpc@23.9'

levante_gcc_netcdf_yac=netcdf-c/4.8.1-openmpi-4.1.2-gcc-11.2.0
levante_gcc_openblas_yac=openblas@0.3.18%gcc@=11.2.0
