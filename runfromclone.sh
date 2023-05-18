#!/bin/bash

### ----- You need to edit these lines to set your ----- ###
### -----  default compiler and python environment ----- ###
module load gcc/11.2.0-gcc-11.2.0
module load python3/2022.01-gcc-11.2.0
source activate /work/mh1126/m300950/pySDenv
### ---------------------------------------------------- ###

### build CLEO (with openMP thread parallelism using Kokkos)
CXX=g++ CC=gcc cmake -S ./ -B ./build -DKokkos_ENABLE_OPENMP=ON -DKokkos_ARCH_NATIVE=ON

### it's a good idea to ensure these directories exist
mkdir ./build/bin
mkdir ./build/share

### generate input files
python ./create_gbxboundariesbinary_script.py ./
python ./create_initsuperdropsbinary_script.py ./ 
python ./create_initthermobinary_script.py ./

### compile and run CLEO
cd build
make clean && make -j 16
./src/runCLEO "../src/config/config.txt" "../libs/claras_SDconstants.hpp"