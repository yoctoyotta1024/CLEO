#!/bin/bash
#SBATCH --job-name=buildgpu
#SBATCH --partition=gpu
#SBATCH --gpus=4
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --mem=30G
#SBATCH --time=00:05:00
#SBATCH --mail-user=clara.bayley@mpimet.mpg.de
#SBATCH --mail-type=FAIL
#SBATCH --account=mh1126
#SBATCH --output=./build/bin/buildgpu_out.%j.out
#SBATCH --error=./build/bin/buildgpu_err.%j.out

### ----- You need to edit these lines to set your ----- ###
### ----- default compiler and python environment   ---- ###
### ----  and paths for CLEO and build directories  ---- ###
module load gcc/11.2.0-gcc-11.2.0
module load nvhpc/23.9-gcc-11.2.0
spack load cmake@3.23.1%gcc
source activate /work/mh1126/m300950/condaenvs/cleoenv
path2CLEO=${HOME}/CLEO/
path2build=$1 # get from command line argument
gxx="g++"
gcc="gcc"
### ---------------------------------------------------- ###

### ------------ choose CUDA compiler ----------- ###
# set path to Kokkos nvcc wrapper (usually Kokkos bin directory of kokkos after installation)
nvcc_wrapper="${HOME}/CLEO/build/_deps/kokkos-src/bin/nvcc_wrapper"

# set nvcc compiler used by Kokkos nvcc wrapper as CUDA_ROOT/bin/nvcc
# NOTE(!) this path should correspond to the loaded nvhpc module. Get path e.g. via 'spack find -p nvhpc@23.9'
CUDA_ROOT="/sw/spack-levante/nvhpc-23.9-xpxqeo/Linux_x86_64/23.9/cuda/"
### ---------------------------------------------------- ###

### ------------ choose Kokkos configuration ----------- ###
# flags for serial kokkos
kokkosflags="-DKokkos_ARCH_NATIVE=ON -DKokkos_ARCH_AMPERE80=ON -DKokkos_ENABLE_SERIAL=ON"

# flags for host parallelism (e.g. using OpenMP)
kokkoshost="-DKokkos_ENABLE_OPENMP=ON"

# flags for device parallelism (e.g. on gpus)
kokkosdevice="-DKokkos_ENABLE_CUDA=ON -DKokkos_ENABLE_CUDA_LAMBDA=O -DKokkos_ENABLE_CUDA_CONSTEXPR=ON \
-DCUDA_ROOT=${CUDA_ROOT} -DKOKKOS_NVCC_WRAPPER=${nvcc_wrapper}"
### ---------------------------------------------------- ###

### ------------ choose extra compiler flags ----------- ###
# flags="-g -O0 -mpc64"                                        # correctness
flags="-O3"                                                    # performance
### ---------------------------------------------------- ###

### ------------ build and compile with cmake ---------- ###
echo "CXX=${gxx} CC=${gcc}"
echo "CUDA=${CUDA_ROOT}/bin/nvcc (via Kokkos nvcc wrapper)"
echo "KOKKOS NVCC WRAPPER=${nvcc_wrapper}"
echo "BUILD_DIR: ${path2build}"
echo "KOKKOS_FLAGS: ${kokkosflags}"
echo "KOKKOS_DEVICE_PARALLELISM: ${kokkosdevice}"
echo "KOKKOS_HOST_PARALLELISM: ${kokkoshost}"
echo "CXX_COMPILER_FLAGS: ${flags}"

# build then compile in parallel
cmake -DCMAKE_CXX_COMPILER=${gxx} \
    -DCMAKE_CC_COMPILER=${gcc} \
    -DCMAKE_CXX_FLAGS="${flags}" \
    -S ${path2CLEO} -B ${path2build} \
    ${kokkosflags} ${kokkosdevice} ${kokkoshost} && \
    cmake --build ${path2build} --parallel

# ensure these directories exist (it's a good idea for later use)
mkdir ${path2build}bin
mkdir ${path2build}share
### ---------------------------------------------------- ###
