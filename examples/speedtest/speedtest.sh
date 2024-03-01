#!/bin/bash
#SBATCH --job-name=speedtest
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --gpus=4
#SBATCH --ntasks-per-node=128
#SBATCH --mem=30G
#SBATCH --time=00:30:00
#SBATCH --mail-user=clara.bayley@mpimet.mpg.de
#SBATCH --mail-type=FAIL
#SBATCH --account=mh1126
#SBATCH --output=./speedtest_out.%j.out
#SBATCH --error=./speedtest_err.%j.out

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
configfile=${path2CLEO}/examples/speedtest/src/config/speedtest_config.txt

python=/work/mh1126/m300950/condaenvs/superdropsenv/bin/python
gxx="g++"
gcc="gcc"
cuda="nvc++"
### ---------------------------------------------------- ###

### ------------ choose Kokkos configuration ----------- ###
kokkosflags="-DKokkos_ARCH_NATIVE=ON -DKokkos_ARCH_AMPERE80=ON -DKokkos_ENABLE_SERIAL=ON"
export OMP_PROC_BIND=spread
export OMP_PLACES=threads

mkdir ${path2build}/bin
### ---------------------------------------------------- ###

### ------------ build gpu CUDA + cpu OpenMP parallelism ----------- ###
buildtype="gpus_cpus"
kokkoshost="-DKokkos_ENABLE_OPENMP=ON"
kokkosdevice="-DKokkos_ENABLE_CUDA=ON -DKokkos_ENABLE_CUDA_LAMBDA=ON"
buildcmd="CXX=${gxx} CC=${gcc} CUDA=${cuda} cmake -S ${path2CLEO} -B " \
  | "${path2build}${buildtype}/ ${kokkosflags} ${kokkoshost} ${kokkosdevice}"
echo ${buildcmd}
CXX=${gxx} CC=${gcc} CUDA=${cuda} cmake -S ${path2CLEO} -B ${path2build}${buildtype}"/" ${kokkosflags} ${kokkoshost} ${kokkosdevice}

### run test for gpus CUDA + cpus OpenMP
mkdir ${path2build}${buildtype}/bin
mkdir ${path2build}${buildtype}/share

${python} speedtest.py ${path2CLEO} ${path2build}${buildtype}"/" ${configfile} ${path2build}"/bin/" ${buildtype}
### ---------------------------------------------------------------- ###

### ---------------- build cpu OpenMP parallelism ------------------ ###
buildtype="cpus"
kokkoshost="-DKokkos_ENABLE_OPENMP=ON"
buildcmd="CXX=${gxx} CC=${gcc} cmake -S ${path2CLEO} -B" \
  | "${path2build}${buildtype}/ ${kokkosflags} ${kokkoshost}"
echo ${buildcmd}
CXX=${gxx} CC=${gcc} cmake -S ${path2CLEO} -B ${path2build}${buildtype}"/" ${kokkosflags} ${kokkoshost}

### run test for cpus OpenMP
mkdir ${path2build}${buildtype}/bin
mkdir ${path2build}${buildtype}/share

${python} speedtest.py ${path2CLEO} ${path2build}${buildtype}"/" ${configfile} ${path2build}"/bin/" ${buildtype}
### ---------------------------------------------------------------- ###

### -------------------------- build serial ------------------------ ###
buildtype="serial"
kokkoshost="-DKokkos_ENABLE_OPENMP=OFF"
buildcmd="CXX=${gxx} CC=${gcc} cmake -S ${path2CLEO} -B" \
  | "${path2build}${buildtype}/ ${kokkosflags} ${kokkoshost}"
echo ${buildcmd}
CXX=${gxx} CC=${gcc} cmake -S ${path2CLEO} -B ${path2build}${buildtype}"/" ${kokkosflags} ${kokkoshost}

### run test in serial
mkdir ${path2build}${buildtype}/bin
mkdir ${path2build}${buildtype}/share

${python} speedtest.py ${path2CLEO} ${path2build}${buildtype}"/" ${configfile} ${path2build}"/bin/" ${buildtype}
### ---------------------------------------------------------------- ###
