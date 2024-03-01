#!/bin/bash
#SBATCH --job-name=rain1d
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --mem=30G
#SBATCH --time=00:10:00
#SBATCH --mail-user=clara.bayley@mpimet.mpg.de
#SBATCH --mail-type=FAIL
#SBATCH --account=mh1126
#SBATCH --output=./rain1d_out.%j.out
#SBATCH --error=./rain1d_err.%j.out

### ----- You need to edit these lines to set your ----- ###
### ----- default compiler and python environment   ---- ###
### ----  and paths for CLEO and build directories  ---- ###
module load gcc/11.2.0-gcc-11.2.0
module load python3/2022.01-gcc-11.2.0
module load nvhpc/23.7-gcc-11.2.0
spack load cmake@3.23.1%gcc
source activate /work/mh1126/m300950/condaenvs/superdropsenv

path2CLEO=${HOME}/CLEO/
path2build=${HOME}/CLEO/build/
configfile=${path2CLEO}/examples/rainshaft1d/src/config/rain1d_config.txt

python=/work/mh1126/m300950/condaenvs/superdropsenv/bin/python
gxx="g++"
gcc="gcc"
cuda="nvc++"
### ---------------------------------------------------- ###

### ------------ choose Kokkos configuration ----------- ###
kokkosflags="-DKokkos_ARCH_NATIVE=ON -DKokkos_ARCH_AMPERE80=ON -DKokkos_ENABLE_SERIAL=ON" # serial kokkos
kokkoshost="-DKokkos_ENABLE_OPENMP=ON"                                          # flags for host parallelism (e.g. using OpenMP)
kokkosdevice="-DKokkos_ENABLE_CUDA=ON -DKokkos_ENABLE_CUDA_LAMBDA=ON"           # flags for device parallelism (e.g. using CUDA)
### ---------------------------------------------------- ###

### ------------------------ build --------------------- ###
### build CLEO using cmake (with openMP thread parallelism through Kokkos)
buildcmd="CXX=${gxx} CC=${gcc} CUDA=${cuda} cmake -S ${path2CLEO} -B ${path2build} ${kokkosflags} ${kokkoshost} ${kokkosdevice}"
echo ${buildcmd}
CXX=${gxx} CC=${gcc} CUDA=${cuda} cmake -S ${path2CLEO} -B ${path2build} ${kokkosflags} ${kokkoshost} ${kokkosdevice}

export OMP_PROC_BIND=spread
export OMP_PLACES=threads

### ensure these directories exist (it's a good idea for later use)
mkdir ${path2build}bin
mkdir ${path2build}share
### ---------------------------------------------------- ###

### ------------------- compile & run ------------------ ###
### generate input files and run 1-D rainshaft example
${python} rainshaft1d.py ${path2CLEO} ${path2build} ${configfile}

### ---------------------------------------------------- ###
