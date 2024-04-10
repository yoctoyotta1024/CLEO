#!/bin/bash
#SBATCH --job-name=build_cleocoupledsdm
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --mem=30G
#SBATCH --time=00:05:00
#SBATCH --mail-user=clara.bayley@mpimet.mpg.de
#SBATCH --mail-type=FAIL
#SBATCH --account=mh1126
#SBATCH --output=./build/bin/build_cleocoupledsdm_out.%j.out
#SBATCH --error=./build/bin/build_cleocoupledsdm_err.%j.out

### ------------------ input parameters ---------------- ###
### ----- You need to edit these lines to specify ------ ###
### ------ (your environment and) directory paths ------ ###
spack load cmake@3.23.1%gcc

buildtype=$1
executable="cleocoupledsdm"
path2CLEO=${HOME}/CLEO/
path2build=$2 # get from command line argument

if [ "${path2build}" == "" ]
then
  path2build=${HOME}/CLEO/build/
fi
### ---------------------------------------------------- ###

### ----------------- build executable --------------- ###
buildcmd="${path2CLEO}/scripts/bash/build_cleo.sh ${buildtype} ${path2CLEO} ${path2build}"
echo ${buildcmd}
${buildcmd}
### ---------------------------------------------------- ###

### ----------------- compile executable --------------- ###
cd ${path2build} && make clean

compilecmd="${path2CLEO}/scripts/bash/compile_cleo.sh ${path2build} ${executable}"
echo ${compilecmd}
${compilecmd}
### ---------------------------------------------------- ###
