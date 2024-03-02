#!/bin/bash
#SBATCH --job-name=buildcpu
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --mem=30G
#SBATCH --time=00:05:00
#SBATCH --mail-user=clara.bayley@mpimet.mpg.de
#SBATCH --mail-type=FAIL
#SBATCH --account=mh1126
#SBATCH --output=./build/bin/buildcpu_out.%j.out
#SBATCH --error=./build/bin/buildcpu_err.%j.out

### ---------------------------------------------------- ###
### ------- You MUST edit these lines to set your ------ ###
### --- default compiler(s) (and python environment) --- ###
### ----  and paths for CLEO and build directories  ---- ###
### ---------------------------------------------------- ###
module load gcc/11.2.0-gcc-11.2.0
spack load cmake@3.23.1%gcc
source activate /work/mh1126/m300950/condaenvs/cleoenv
path2CLEO=${HOME}/CLEO/
path2build=$1             # get from command line argument(s)
gxx="g++"
gcc="gcc"
### ---------------------------------------------------- ###

### ---------------------------------------------------- ###
### ------- You can optionally edit the following ------ ###
### -------- lines to customise your compiler(s) ------- ###
###  ------------ and build configuration  ------------- ###
### ---------------------------------------------------- ###

### ------------ choose extra compiler flags ----------- ###
# flags="-g -O0 -mpc64"                                        # correctness
flags="-O3"                                                    # performance
### ---------------------------------------------------- ###

### ------------ choose Kokkos configuration ----------- ###
# flags for serial kokkos
kokkosflags="-DKokkos_ARCH_NATIVE=ON -DKokkos_ARCH_AMPERE80=ON -DKokkos_ENABLE_SERIAL=ON"

# flags for host parallelism (e.g. using OpenMP)
kokkoshost="-DKokkos_ENABLE_OPENMP=ON"

# flags for device parallelism (e.g. on gpus)
kokkosdevice=""
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
