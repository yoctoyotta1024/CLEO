#!/bin/bash

### -------------- GCC compiler(s) Packages ------------ ###
juwels_gcc=GCC/13.3.0 # module load
juwels_gcc_cmake=CMake/3.29.3 # module load
juwels_gcc_mpi=ParaStationMPI/5.10.0-1 # module load
juwels_gxx_compiler="/p/software/default/stages/2025/software/psmpi/5.10.0-1-GCC-13.3.0/bin/mpic++"
juwels_gcc_compiler="/p/software/default/stages/2025/software/psmpi/5.10.0-1-GCC-13.3.0/bin/mpicc"
juwels_gcc_gsl=GSL/2.8 # module load (for MPTRAC)
juwels_gcc_netcdf=netCDF/4.9.2 # module load (for MPTRAC)
### ---------------------------------------------------- ###

### ------------- Intel compiler(s) Packages ------------ ###
juwels_intel=Intel/2024.2.0 # module load
juwels_intel_cmake=CMake/3.29.3 # module load
juwels_intel_mpi=ParaStationMPI/5.10.0-1 # module load
juwels_icpc_compiler="/p/software/default/stages/2025/software/psmpi/5.10.0-1-intel-compilers-2024.2.0/bin/mpic++"
juwels_icc_compiler="/p/software/default/stages/2025/software/psmpi/5.10.0-1-intel-compilers-2024.2.0/bin/mpicc"
### ---------------------------------------------------- ###
