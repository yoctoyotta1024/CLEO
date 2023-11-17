#!/bin/bash
#SBATCH --job-name=adiabaticparcel
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --mem=30G
#SBATCH --time=00:30:00
#SBATCH --mail-user=clara.bayley@mpimet.mpg.de
#SBATCH --mail-type=FAIL
#SBATCH --account=mh1126
#SBATCH --output=./build/bin/adiabaticparcel_out.%j.out
#SBATCH --error=./build/bin/adiabaticparcel_err.%j.out

### ----- You need to edit these lines to set your ----- ###
### ----- default compiler and python environment   ---- ###
### ----  and paths for CLEO and build directories  ---- ###
module load gcc/11.2.0-gcc-11.2.0
source activate /work/mh1126/m300950/condaenvs/cleoenv
path2CLEO=${HOME}/CLEO/
path2build=${HOME}/CLEO/build/
configfile=${HOME}/examples/adiabaticparcel/src/config/config.txt 
python=python
gxx="g++"
gcc="gcc"

### build CLEO using cmake (with openMP thread parallelism through Kokkos)
kokkosflags="-DKokkos_ARCH_NATIVE=ON -DKokkos_ENABLE_SERIAL=ON -DKokkos_ENABLE_OPENMP=ON"  # openMP parallelism enabled
CXX=${gxx} CC=${gcc} cmake -S ${path2CLEO} -B ${path2build} ${kokkosflags}

### ensure these directories exist (it's a good idea for later use)
mkdir ${path2build}bin
mkdir ${path2build}share

### generate input files and run adiabatic parcel example
${python} condtest.py ${path2CLEO} ${path2build} ${configfile}