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
### -----  default compiler and python environment ----- ###
module load gcc/11.2.0-gcc-11.2.0
module load python3/2022.01-gcc-11.2.0
source activate /work/mh1126/m300950/superdropsV2
### ---------------------------------------------------- ###

### build CLEO (with openMP thread parallelism using Kokkos)
CXX=g++ CC=gcc cmake -S ./ -B ./build -DKokkos_ENABLE_OPENMP=ON -DKokkos_ARCH_NATIVE=ON

### it's a good idea to ensure these directories exist
mkdir ./build/bin
mkdir ./build/share

### generate input files
# path2CLEO = "/home/m/m300950/CLEO/"
# path2CLEO = "/Users/yoctoyotta1024/Documents/b1_springsummer2023/CLEO/"
# path2build = "/home/m/m300950/CLEO/build/"
# path2build = "/Users/yoctoyotta1024/Documents/b1_springsummer2023/CLEO/build/"
python ./create_gbxboundariesbinary_script.py ./
python ./create_initsuperdropsbinary_script.py ./ 
python ./create_initthermobinary_script.py ./

### compile and run CLEO
cd build
make clean && make -j 16
./src/runCLEO "../src/config/config.txt" "../libs/claras_SDconstants.hpp"