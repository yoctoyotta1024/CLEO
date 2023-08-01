#!/bin/bash
#SBATCH --job-name=quickrun
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --mem=30G
#SBATCH --time=00:30:00
#SBATCH --mail-user=clara.bayley@mpimet.mpg.de
#SBATCH --mail-type=FAIL
#SBATCH --account=mh1126
#SBATCH --output=./build/bin/quickrun_out.%j.out
#SBATCH --error=./build/bin/quickrun_err.%j.out

### ----- You need to edit these lines to set your ----- ###
### ----- default compiler and python environment   ---- ###
### ----  and paths for CLEO and build directories  ---- ###
module load gcc/11.2.0-gcc-11.2.0
module load python3/2022.01-gcc-11.2.0
source activate /work/mh1126/m300950/condaenvs/cleoenv 
path2CLEO=${HOME}/CLEO/
path2build=${HOME}/CLEO/example/build/
configfile=${HOME}/CLEO/example/exmpl_config.txt
python=python
gxx="g++"
gcc="gcc"

# path2CLEO=${HOME}/Documents/b1_springsummer2023/CLEO/
# #path2build=${HOME}/CLEO/example/build/                            ### TODO: correct path for my home dir
# #configfile=${HOME}/CLEO/example/exampleconfig.txt                ### TODO: correct path for my home dir
# python=${HOME}/opt/anaconda3/envs/superdropsV2/bin/python
# gxx="g++-13"
# gcc="gcc-13"
## ---------------------------------------------------- ###

### build CLEO using cmake (with openMP thread parallelism through Kokkos)
kokkosflags="-DKokkos_ARCH_NATIVE=ON -DKokkos_ENABLE_SERIAL=ON -DKokkos_ENABLE_OPENMP=ON"  # openMP parallelism enabled
CXX=${gxx} CC=${gcc} cmake -S ${path2CLEO} -B ${path2build} ${kokkosflags}

### ensure these directories exist (it's a good idea for later use)
mkdir ${path2build}bin
mkdir ${path2build}share

### generate input files
${python} exmpl_createinputbinaries.py ${path2CLEO} ${path2build} ${configfile}

### compile CLEO
cd ${path2build} && pwd
make clean && make -j 16

### run CLEO
export OMP_PROC_BIND=spread
export OMP_PLACES=threads
runcmd="${path2build}/src/runCLEO ${configfile} ${path2CLEO}libs/claras_SDconstants.hpp"
echo ${runcmd}
${runcmd}