#!/bin/bash
#SBATCH --job-name=compile_cleo
#SBATCH --partition=gpu
#SBATCH --gpus=4
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --mem=30G
#SBATCH --time=00:30:00
#SBATCH --mail-user=clara.bayley@mpimet.mpg.de
#SBATCH --mail-type=FAIL
#SBATCH --account=mh1126
#SBATCH --output=./compile_cleo_out.%j.out
#SBATCH --error=./compile_cleo_err.%j.out

### ------------- PLEASE NOTE: this script assumes you ------------- ###
### ----------- have already built CLEO in "path2build" ------------ ###
### -------------------  directory using cmake  -------------------- ###

### ------------------ input parameters ---------------- ###
### ----- You need to edit these lines to specify ------ ###
### ----- (your environment and) build directory ------- ###
### -------- path and executable to compile ------------ ###
spack load cmake@3.23.1%gcc

path2build=$1   # get from command line argument
executable=$2   # get from command line argument
### ---------------------------------------------------- ###

### ----------------- compile executable --------------- ###
echo "path to build directory: ${path2build}"
echo "executable: ${executable}"

cd ${path2build}
make clean
make -j 128 ${executable}
### ---------------------------------------------------- ###
