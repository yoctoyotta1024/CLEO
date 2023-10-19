#!/bin/bash
#SBATCH --job-name=gpubuildCLEO
#SBATCH --partition=gpu
#SBATCH --gpus=4
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
module load python3/2022.01-gcc-11.2.0
module load nvhpc/23.7-gcc-11.2.0
spack load cmake@3.23.1%gcc
source activate /work/mh1126/m300950/condaenvs/cleoenv 
path2CLEO=${HOME}/testCLEOfire/
path2build=${HOME}/testCLEOfire/build/
python=python
gxx="nvc++"
gcc="nvcc"
### ---------------------------------------------------- ###

### ------------------ build_compile.sh ---------------- ###
### add c++ library used by nvhpc compiler to runtime path
export LD_LIBRARY_PATH=/sw/spack-levante/gcc-11.2.0-bcn7mb/lib64

### build CLEO using cmake (with optional thread parallelism through Kokkos)
kokkosflags="-DKokkos_ARCH_NATIVE=ON -DKokkos_ARCH_AMPERE80=ON -DKokkos_ENABLE_SERIAL=ON -DKokkos_ENABLE_CUDA=ON -DKokkos_ENABLE_CUDA_LAMBDA=ON"  # CUDA parallelism enabled
CXX=${gxx} CC=${gcc} cmake -S ${path2CLEO} -B ${path2build} ${kokkosflags}

### compile CLEO
cd ${path2build} && pwd 
make clean && make -j 16
### ---------------------------------------------------- ###