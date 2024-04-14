#!/bin/bash
#SBATCH --job-name=cuspbifurc
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --mem=30G
#SBATCH --time=00:05:00
#SBATCH --mail-user=clara.bayley@mpimet.mpg.de
#SBATCH --mail-type=FAIL
#SBATCH --account=mh1126
#SBATCH --output=./cuspbifurc_out.%j.out
#SBATCH --error=./cuspbifurc_err.%j.out

### ---------------------------------------------------- ###
### ------------------ Input Parameters ---------------- ###
### ------ You MUST edit these lines to set your ------- ###
### ---- build type, directories, the executable(s) ---- ###
### -------- to compile, and your python script -------- ###
### ---------------------------------------------------- ###
buildtype="serial"
path2CLEO=${HOME}/CLEO/
path2build=${HOME}/CLEO/build_adia0D/
executables="adia0D"

pythonscript=${path2CLEO}/examples/adiabaticparcel/cuspbifurc.py
configfile=${path2CLEO}/examples/adiabaticparcel/src/config/cuspbifurc_config.txt
### ---------------------------------------------------- ###
### ---------------------------------------------------- ###
### ---------------------------------------------------- ###

### ---------- build, compile and run example ---------- ###
runcmd="${path2CLEO}/examples/run_example.sh ${buildtype} \
  ${path2CLEO} ${path2build} ${executables} ${pythonscript} ${configfile}"
echo ${runcmd}
${runcmd}
### ---------------------------------------------------- ###
