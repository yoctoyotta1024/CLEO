#!/bin/bash
#SBATCH --job-name=as2017
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --mem=30G
#SBATCH --time=00:30:00
#SBATCH --mail-user=clara.bayley@mpimet.mpg.de
#SBATCH --mail-type=FAIL
#SBATCH --account=mh1126
#SBATCH --output=./as2017_out.%j.out
#SBATCH --error=./as2017_err.%j.out

### ----- You need to edit these lines to set your ----- ###
### ----- default compiler and python environment   ---- ###
### ----  and paths for CLEO and build directories  ---- ###
module load gcc/11.2.0-gcc-11.2.0
module load python3/2022.01-gcc-11.2.0
source activate /work/mh1126/m300950/condaenvs/superdropsenv
path2CLEO=${HOME}/CLEO/
path2build=${HOME}/CLEO/build/
configfile=${path2CLEO}/examples/adiabaticparcel/src/config/as2017_config.txt 

python=/work/mh1126/m300950/condaenvs/superdropsenv/bin/python
gxx="g++"
gcc="gcc"
### ---------------------------------------------------- ###

### ------------ choose Kokkos configuration ----------- ###
kokkosflags="-DKokkos_ARCH_NATIVE=ON -DKokkos_ENABLE_SERIAL=ON" # serial kokkos
kokkoshost="-DKokkos_ENABLE_OPENMP=ON"                          # flags for host parallelism (e.g. using OpenMP)
# kokkosdevice=""
kokkosdevice="-DKokkos_ENABLE_CUDA=ON -DKokkos_ENABLE_CUDA_LAMBDA=ON"   # flags for device parallelism (e.g. using CUDA)
### ---------------------------------------------------- ###

### ------------------------ build --------------------- ###
### build CLEO using cmake (with openMP thread parallelism through Kokkos)
buildcmd="CXX=${gxx} CC=${gcc} CUDA=${cuda} cmake -S ${path2CLEO} -B ${path2build} ${kokkosflags} ${kokkoshost}"
echo ${buildcmd}
CXX=${gxx} CC=${gcc} CUDA=${cuda} cmake -S ${path2CLEO} -B ${path2build} ${kokkosflags} ${kokkoshost}

export OMP_PROC_BIND=spread
export OMP_PLACES=threads

### ensure these directories exist (it's a good idea for later use)
mkdir ${path2build}bin
mkdir ${path2build}share
### ---------------------------------------------------- ###

### ------------------- compile & run ------------------ ###
### generate input files and run adiabatic parcel example
${python} as2017.py ${path2CLEO} ${path2build} ${configfile}

### ---------------------------------------------------- ###