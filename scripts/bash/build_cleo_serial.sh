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

path2CLEO=$1    # get from command line argument
path2build=$2   # get from command line argument
### ------------------------------------------------------------------------ ###

### ---------------------------------------------------- ###
### ------- You can optionally edit the following ------ ###
### -------- lines to customise your compiler(s) ------- ###
###  ------------ and build configuration  ------------- ###
### ---------------------------------------------------- ###

### --------- choose C/C++ compiler and flags ---------- ###
CC=${gcc}               # C
CXX=${gxx}              # C++

# CMAKE_CXX_FLAGS="-Werror -Wall -Wextra -pedantic -g -gdwarf-4 -O0 -mpc64"      # correctness and debugging (note -gdwarf-4 not possible for nvc++)
CMAKE_CXX_FLAGS="-Werror -Wall -pedantic -O3"                            # performance
### ---------------------------------------------------- ###

### ------------ choose Kokkos configuration ----------- ###
# flags for serial kokkos
kokkosflags="-DKokkos_ARCH_NATIVE=ON -DKokkos_ENABLE_SERIAL=ON" # serial kokkos

# flags for host parallelism (e.g. using OpenMP)
kokkoshost=""

# flags for device parallelism (e.g. on gpus)
kokkosdevice=""
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
    -DCMAKE_CC_COMPILER=${CC} \
    -DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS}" \
    -S ${path2CLEO} -B ${path2build} \
    ${kokkosflags} ${kokkosdevice} ${kokkoshost}

# ensure these directories exist (it's a good idea for later use)
mkdir -p ${path2build}/bin
mkdir -p ${path2build}/share
### ---------------------------------------------------- ###
