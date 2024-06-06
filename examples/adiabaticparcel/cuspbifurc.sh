#!/bin/bash
#SBATCH --job-name=cuspbifurc
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --mem=30G
#SBATCH --time=00:10:00
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
enableyac=false
executables="adia0D"

pythonscript=${path2CLEO}/examples/adiabaticparcel/cuspbifurc.py
configfile=${path2CLEO}/examples/adiabaticparcel/src/config/cuspbifurc_config.yaml
script_args="${configfile}"
### ---------------------------------------------------- ###
### ---------------------------------------------------- ###
### ---------------------------------------------------- ###

### ---------- build, compile and run example ---------- ###
${path2CLEO}/examples/run_example.sh \
  ${buildtype} ${path2CLEO} ${path2build} ${enableyac} \
  "${executables}" ${pythonscript} "${script_args}"
### ---------------------------------------------------- ###
