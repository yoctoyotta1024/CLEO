#!/bin/bash
#SBATCH --job-name=as2017
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --mem=30G
#SBATCH --time=00:10:00
#SBATCH --mail-user=clara.bayley@mpimet.mpg.de
#SBATCH --mail-type=FAIL
#SBATCH --account=mh1126
#SBATCH --output=./as2017_out.%j.out
#SBATCH --error=./as2017_err.%j.out

### ------ Generic script to build CLEO, compile  ------ ###
### ----- some of its executables and run a python ----- ###
### ------  script e.g. for example(s) on Levante. ----- ###

### ---------------------------------------------------- ###
### ------------------ Input Parameters ---------------- ###
### ------ You MUST edit these lines to set your ------- ###
### ---- build type, directories, the executable to ---- ###
### ----------- compile and your environment ----------- ###
### ---------------------------------------------------- ###
spack load cmake@3.23.1%gcc
module load python3/2022.01-gcc-11.2.0
source activate /work/mh1126/m300950/condaenvs/superdropsenv
python=/work/mh1126/m300950/condaenvs/superdropsenv/bin/python

buildtype="openmp"
path2CLEO=${HOME}/CLEO/
path2build=${HOME}/CLEO/build0/
executable="adia0D"
configfile=${path2CLEO}/examples/adiabaticparcel/src/config/as2017_config.txt
### ---------------------------------------------------- ###
### ---------------------------------------------------- ###
### ---------------------------------------------------- ###

### ---------------------- build CLEO ------------------ ###
buildcmd="${path2CLEO}/scripts/bash/build_cleo.sh ${buildtype} ${path2CLEO} ${path2build}"
echo ${buildcmd}
${buildcmd}
### ---------------------------------------------------- ###

### --------- compile executable from scratch ---------- ###
cd ${path2build} && make clean

compilecmd="${path2CLEO}/scripts/bash/compile_cleo.sh ${buildtype} ${path2build} ${executable}"
echo ${compilecmd}
${compilecmd}
### ---------------------------------------------------- ###

### --------- run model through Python script ---------- ###
export OMP_PROC_BIND=spread
export OMP_PLACES=threads

${python} ${path2CLEO}/examples/adiabaticparcel/as2017.py \
  ${path2CLEO} ${path2build} ${configfile}
### ---------------------------------------------------- ###
