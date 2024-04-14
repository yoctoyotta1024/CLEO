#!/bin/bash
#SBATCH --job-name=as2017
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --mem=30G
#SBATCH --time=00:05:00
#SBATCH --mail-user=clara.bayley@mpimet.mpg.de
#SBATCH --mail-type=FAIL
#SBATCH --account=mh1126
#SBATCH --output=./as2017_out.%j.out
#SBATCH --error=./as2017_err.%j.out

### ---------------------------------------------------- ###
### ------------------ Input Parameters ---------------- ###
### ------ You MUST edit these lines to set your ------- ###
### ---- build type, directories, the executable to ---- ###
### ---------- compile, and your python script --------- ###
### ---------------------------------------------------- ###
buildtype="openmp"
path2CLEO=${HOME}/CLEO/
path2build=${HOME}/CLEO/build_adia0D/
executable="adia0D"

pythonscript=${path2CLEO}/examples/adiabaticparcel/as2017.py
configfile=${path2CLEO}/examples/adiabaticparcel/src/config/as2017_config.txt
### ---------------------------------------------------- ###
### ---------------------------------------------------- ###
### ---------------------------------------------------- ###

### ---------- build, compile and run example ---------- ###
runcmd="${path2CLEO}/examples/run_example.sh ${buildtype} \
  ${path2CLEO} ${path2build} ${executable} ${pythonscript} ${configfile}"
echo ${runcmd}
${runcmd}
### ---------------------------------------------------- ###
