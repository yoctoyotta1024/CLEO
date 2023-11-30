#!/bin/bash

### ----- You need to edit these lines to set your ----- ###
### ----- default compiler and python environment   ---- ###
### ----  and paths for CLEO and build directories  ---- ###
module load gcc/11.2.0-gcc-11.2.0
module load python3/2022.01-gcc-11.2.0
module load nvhpc/23.7-gcc-11.2.0
spack load cmake@3.23.1%gcc
source activate /work/mh1126/m300950/condaenvs/cleoenv 
path2CLEO=${HOME}/CLEO/
configfile=${HOME}/CLEO/src/config/config.txt

python=python
gxx="g++"
gcc="gcc"
cuda="nvc++"
### ---------------------------------------------------- ###

export OMP_PROC_BIND=spread
export OMP_PLACES=threads

cd ${HOME}/CLEO/buildserial/ && pwd 
make -j 64

cd ${HOME}/CLEO/build/ && pwd 
make -j 64

cd ${HOME}/CLEO/buildgpu/ && pwd 
make -j 64