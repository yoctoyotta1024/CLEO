#!/bin/bash
#SBATCH --job-name=roughCLEO
#SBATCH --partition=gpu
#SBATCH --gpus=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --mem=30G
#SBATCH --time=00:05:00
#SBATCH --mail-user=clara.bayley@mpimet.mpg.de
#SBATCH --mail-type=FAIL
#SBATCH --account=mh1126
#SBATCH --output=./gpubuildCLEO_out.%j.out
#SBATCH --error=./gpubuildCLEO_err.%j.out

### ----- You need to edit these lines to set your ----- ###
### ----- default compiler and python environment   ---- ###
### ----  and paths for CLEO and build directories  ---- ###
module load gcc/11.2.0-gcc-11.2.0
module load nvhpc/23.7-gcc-11.2.0
spack load cmake@3.23.1%gcc
source activate /work/mh1126/m300950/condaenvs/cleoenv
path2build=${HOME}/CLEO/roughpaper/build/
gxx="g++"
gcc="gcc"
cuda="nvcc"
### ---------------------------------------------------- ###

### ------------ choose Kokkos configuration ----------- ###

kokkosflags="-DKokkos_ARCH_NATIVE=ON -DKokkos_ARCH_AMPERE80=ON -DKokkos_ENABLE_SERIAL=ON"                 # serial kokkos
# kokkosdevice="-DKokkos_ENABLE_CUDA=ON -DKokkos_ENABLE_CUDA_LAMBDA=ON -DKokkos_ENABLE_CUDA_CONSTEXPR=ON" # flags for device parallelism (e.g. on gpus)
# kokkoshost="-DKokkos_ENABLE_OPENMP=ON"                                                                  # flags for host parallelism (e.g. using OpenMP)
### ---------------------------------------------------- ###

### ------------ choose extra compiler flags ----------- ###
flags="-g -O0 -mpc64"                                            # correctness
# flags="-O3"                                                    # performance
### ---------------------------------------------------- ###

### ------------------ build_compile.sh ---------------- ###
### build CLEO using cmake (with optional thread parallelism through Kokkos)
echo "CXX=${gxx} CC=${gcc} CUDA=${cuda}"
echo "BUILD_DIR: ${path2build}"
echo "KOKKOS_FLAGS: ${kokkosflags} ${kokkosdevice} ${kokkoshost}"
echo "CXX_COMPILER_FLAGS: ${flags}"

CXX=${gxx} CC=${gcc} CUDA=${cuda} cmake -DCMAKE_CXX_COMPILER=${gxx} \
                                      -DCMAKE_CC_COMPILER=${gcc} \
                                      -DCMAKE_CUDA_COMPILER=${cuda} \
                                      -DCMAKE_CXX_FLAGS="${flags}" \
                                      -S ../ -B ${path2build} \
                                      ${kokkosflags} ${kokkosdevice} ${kokkoshost} && \
                                      cmake --build ${path2build}  --parallel

### compile CLEO
cd ${path2build} && pwd
make -j 16
### ---------------------------------------------------- ###
