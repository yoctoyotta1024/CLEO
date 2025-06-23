#!/bin/bash

### -------------- GCC compiler(s) Packages ------------ ###
levante_gcc="gcc/11.2.0-gcc-11.2.0" # bcn7mbu # module load
levante_gcc_cmake="cmake@3.26.3%gcc@=11.2.0/fuvwuhz" # spack load
levante_gcc_openmpi="openmpi@4.1.2%gcc@11.2.0" # spack load
levante_gxx_compiler="/sw/spack-levante/openmpi-4.1.2-mnmady/bin/mpic++"
levante_gcc_compiler="/sw/spack-levante/openmpi-4.1.2-mnmady/bin/mpicc"
levante_gcc_cuda="cuda@12.2.0%gcc@=11.2.0"  # spack load
levante_gcc_cuda_root="/sw/spack-levante/cuda-12.2.0-2ttufp/" # [cuda_root]/bin/nvcc") (can get hint for correct path via 'spack find -p nvhpc@23.9')

### specific gcc compiler compatible packages for YAC installation and usage
levante_gcc_netcdf_yac="netcdf-c/4.8.1-openmpi-4.1.2-gcc-11.2.0" # module load
levante_gcc_openblas_yac="openblas@0.3.18%gcc@=11.2.0" # spack load
levante_gcc_fyaml_root="/sw/spack-levante/libfyaml-0.7.12-fvbhgo" # match `spack location -i libfyaml`
levante_gcc_fyamllib="${levante_gcc_fyaml_root}/lib"
### specific packages for YAC installation only
levante_gcc_f90_compiler="/sw/spack-levante/openmpi-4.1.2-mnmady/bin/mpifort"
levante_gcc_netcdf_root="/sw/spack-levante/netcdf-c-4.8.1-6qheqr" # match `module show ${netcdf}`
### ---------------------------------------------------- ###

### ------------- Intel compiler(s) Packages ------------ ###
levante_intel="intel-oneapi-mpi/2021.5.0-gcc-11.2.0" # module load
levante_intel_cmake="cmake@3.23.1%oneapi" # spack load
levante_intel_openmpi="openmpi@4.1.2%intel@=2021.5.0/yfwe6t4" # spack load
levante_icpc_compiler="/sw/spack-levante/openmpi-4.1.2-yfwe6t/bin/mpic++"
levante_icc_compiler="/sw/spack-levante/openmpi-4.1.2-yfwe6t/bin/mpicc"

### specific gcc compiler compatible packages for YAC installation and usage
levante_intel_netcdf_yac="netcdf-c/4.8.1-intel-oneapi-mpi-2021.5.0-gcc-11.2.0" # module load
levante_intel_openblas_yac="openblas@0.3.18%intel@=2021.5.0" # spack load
levante_intel_fyaml_root="/sw/spack-levante/libfyaml-0.7.12-fvbhgo" # match `spack location -i libfyaml`
levante_intel_fyamllib="${levante_intel_fyaml_root}/lib"
### specific packages for YAC installation only
levante_intel_ifort_compiler="/sw/spack-levante/openmpi-4.1.2-yfwe6t/bin/mpifort"
levante_intel_netcdf_root="/sw/spack-levante/netcdf-c-4.8.1-xguss5" # match `module show ${netcdf}`
### ---------------------------------------------------- ###
