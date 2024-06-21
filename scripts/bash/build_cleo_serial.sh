#!/bin/bash
#SBATCH --job-name=serialbuild
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --mem=30G
#SBATCH --time=00:05:00
#SBATCH --mail-user=clara.bayley@mpimet.mpg.de
#SBATCH --mail-type=FAIL
#SBATCH --account=mh1126
#SBATCH --output=./build/bin/serialbuild_out.%j.out
#SBATCH --error=./build/bin/serialbuild_err.%j.out

### ------------------------------------------------------------------------ ###
### ------- You MUST edit these lines to set your default compiler(s) ------ ###
### --------- and optionally your environment, path to CLEO and the -------- ###
### ----------------------- desired build directory  ----------------------- ###
### ------------------------------------------------------------------------ ###
module load gcc/11.2.0-gcc-11.2.0
spack load cmake@3.23.1%gcc
gxx="/sw/spack-levante/gcc-11.2.0-bcn7mb/bin/g++"
gcc="/sw/spack-levante/gcc-11.2.0-bcn7mb/bin/gcc"

path2CLEO=$1    # required
path2build=$2   # required
enableyac=$3    # required
yacyaxtroot=$4      # required if enableyac=true

yac_openmpi=openmpi/4.1.2-gcc-11.2.0 # required if enableyac=true (must match gcc)
yac_netcdf=netcdf-c/4.8.1-openmpi-4.1.2-gcc-11.2.0 # required if enableyac=true (must match gcc & openmp)
yac_openblas=openblas@0.3.18%gcc@=11.2.0 # required if enableyac=true (must match gcc)
### ------------------------------------------------------------------------ ###

### ---------------------------------------------------- ###
### ------- You can optionally edit the following ------ ###
### -------- lines to customise your compiler(s) ------- ###
###  ------------ and build configuration  ------------- ###
### ---------------------------------------------------- ###

### --------- choose C/C++ compiler and flags ---------- ###
CC=${gcc}               # C
CXX=${gxx}              # C++

## for correctness and debugging (note -gdwarf-4 not possible for nvc++) use:
# CMAKE_CXX_FLAGS="-Werror -Wno-unused-parameter -Wall -Wextra -pedantic -g -gdwarf-4 -O0 -mpc64"
# for performance use:
CMAKE_CXX_FLAGS="-Werror -Wall -pedantic -O3"
### ---------------------------------------------------- ###

### ------------ choose Kokkos configuration ----------- ###
# flags for serial kokkos
kokkosflags="-DKokkos_ARCH_NATIVE=ON -DKokkos_ENABLE_SERIAL=ON" # serial kokkos

# flags for host parallelism (e.g. using OpenMP)
kokkoshost=""

# flags for device parallelism (e.g. on gpus)
kokkosdevice=""
### ---------------------------------------------------- ###

### ------------------ choose YAC build ---------------- ###
if ! [ "${enableyac}" == "true" ]
then
    yacflags="-DENABLE_YAC_COUPLING=OFF"

else
    module load ${yac_openmpi} ${yac_netcdf}
    spack load ${yac_openblas}
    yacflags="-DENABLE_YAC_COUPLING=ON -DYAXT_ROOT=${yacyaxtroot}/yaxt -DYAC_ROOT=${yacyaxtroot}/yac"
    yacmodule="${path2CLEO}/libs/coupldyn_yac/cmake"
    echo "YAC FLAGS: ${yacflags} -DCMAKE_MODULE_PATH=${yacmodule}"
fi
### ---------------------------------------------------- ###

### ---------------- build CLEO with cmake ------------- ###
echo "CXX_COMPILER=${CXX} CC_COMPILER=${CC}"
echo "CLEO_DIR: ${path2CLEO}"
echo "BUILD_DIR: ${path2build}"
echo "KOKKOS_FLAGS: ${kokkosflags}"
echo "KOKKOS_DEVICE_PARALLELISM: ${kokkosdevice}"
echo "KOKKOS_HOST_PARALLELISM: ${kokkoshost}"
echo "CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}"

cmake -DCMAKE_CXX_COMPILER=${CXX} \
    -DCMAKE_C_COMPILER=${CC} \
    -DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS}" \
    -DCMAKE_MODULE_PATH=${yacmodule} \
    -S ${path2CLEO} -B ${path2build} \
    ${kokkosflags} ${kokkosdevice} ${kokkoshost} \
    ${yacflags}

# ensure these directories exist (it's a good idea for later use)
mkdir -p ${path2build}/bin
mkdir -p ${path2build}/share
### ---------------------------------------------------- ###
