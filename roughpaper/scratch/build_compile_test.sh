#!/bin/bash
#SBATCH --job-name=roughCLEO
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --gpus=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=10G
#SBATCH --time=00:05:00
#SBATCH --mail-user=clara.bayley@mpimet.mpg.de
#SBATCH --mail-type=FAIL
#SBATCH --account=bm1183
#SBATCH --output=./roughCLEO_out.%j.out
#SBATCH --error=./roughCLEO_err.%j.out

### ---------------------------------------------------- ###
### ------- You MUST edit these lines to set your ------ ###
### --- default compiler(s) (and python environment) --- ###
### ----  and paths for CLEO and build directories  ---- ###
### ---------------------------------------------------- ###
dobuild=$1 # == "build" or otherwise

source /etc/profile
module purge
spack unload --all
module load gcc/11.2.0-gcc-11.2.0 openmpi/4.1.2-gcc-11.2.0
# module load gcc/11.2.0-gcc-11.2.0 openmpi/4.1.2-gcc-11.2.0 nvhpc/23.9-gcc-11.2.0 # (CUDA)
spack load cmake@3.23.1%gcc

gxx="/sw/spack-levante/openmpi-4.1.2-mnmady/bin/mpic++"
gcc="/sw/spack-levante/openmpi-4.1.2-mnmady/bin/mpicc"

path2CLEO=${HOME}/CLEO/
path2build=${HOME}/CLEO/roughpaper/scratch/build/
### ---------------------------------------------------- ###

### ---------------------------------------------------- ###
### ------- You can optionally edit the following ------ ###
### -------- lines to customise your compiler(s) ------- ###
###  ------------ and build configuration  ------------- ###
### ---------------------------------------------------- ###

### --------- choose C/C++ compiler and flags ---------- ###
CC=${gcc}               # C
CXX=${gxx}              # C++

## for correctness and debugging (note -gdwarf-4 not possible for nvc++) use:
CMAKE_CXX_FLAGS="-Werror -Wno-unused-parameter -Wall -Wextra -pedantic -g -gdwarf-4 -O0 -mpc64"
# for performance use:
# CMAKE_CXX_FLAGS="-Werror -Wall -pedantic -O3"
### ---------------------------------------------------- ###

### --------------- choose CUDA compiler --------------- ###
# set nvcc compiler used by Kokkos nvcc wrapper as CLEO_CUDA_ROOT/bin/nvcc
# NOTE(!) this path should correspond to the loaded nvhpc module.
# you can get a clue for the correct path e.g. via 'spack find -p nvhpc@23.9'
CUDA_ROOT="/sw/spack-levante/nvhpc-23.9-xpxqeo/Linux_x86_64/23.9/cuda/"

# set path to Kokkos nvcc wrapper (usually Kokkos bin directory of kokkos after installation)
NVCC_WRAPPER_DEFAULT_COMPILER=${gxx}
### ---------------------------------------------------- ###

### ------------ choose Kokkos configuration ----------- ###
# flags for serial kokkos
kokkosflags="-DKokkos_ARCH_NATIVE=ON -DKokkos_ARCH_AMPERE80=ON -DKokkos_ENABLE_SERIAL=ON"

# flags for host parallelism (e.g. using OpenMP)
kokkoshost="-DKokkos_ENABLE_OPENMP=ON"

# # flags for device parallelism (e.g. on gpus)
# kokkosdevice="-DKokkos_ENABLE_CUDA=ON \
# -DKokkos_ENABLE_CUDA_CONSTEXPR=ON -DKokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE=ON \
# -DCUDA_ROOT=${CUDA_ROOT} -DNVCC_WRAPPER_DEFAULT_COMPILER=${NVCC_WRAPPER_DEFAULT_COMPILER}"
### ---------------------------------------------------- ###

### ------------ build and compile with cmake ---------- ###
echo "CXX_COMPILER=${CXX} CC_COMPILER=${CC}"
echo "CUDA=${CUDA_ROOT}/bin/nvcc (via Kokkos nvcc wrapper)"
echo "NVCC_WRAPPER_DEFAULT_COMPILER=${NVCC_WRAPPER_DEFAULT_COMPILER}"
echo "BUILD_DIR: ${path2build}"
echo "KOKKOS_FLAGS: ${kokkosflags}"
echo "KOKKOS_DEVICE_PARALLELISM: ${kokkosdevice}"
echo "KOKKOS_HOST_PARALLELISM: ${kokkoshost}"
echo "CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}"

# delete any existing "test" executable
rm ${path2build}/roughpaper/scratch/test

# build then compile in parallel
if [[ ${dobuild} == "build" ]];
then
  cmake -DCMAKE_CXX_COMPILER=${CXX} \
      -DCMAKE_C_COMPILER=${CC} \
      -DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS}" \
      -S ${path2CLEO} -B ${path2build} \
      ${kokkosflags} ${kokkosdevice} ${kokkoshost} && \
      cmake --build ${path2build} --target test --parallel
else
  cd ${path2build} && make -j 64 test
fi

# good settings for Kokkos OpenMP at runtime
export OMP_PROC_BIND=spread
export OMP_PLACES=threads

${path2CLEO}/roughpaper/scratch/build/roughpaper/scratch/test
### ---------------------------------------------------- ###
