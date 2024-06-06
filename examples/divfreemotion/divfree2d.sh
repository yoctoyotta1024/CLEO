#!/bin/bash
#SBATCH --job-name=divfree2d
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --gpus=4
#SBATCH --ntasks-per-node=128
#SBATCH --mem=30G
#SBATCH --time=00:05:00
#SBATCH --mail-user=clara.bayley@mpimet.mpg.de
#SBATCH --mail-type=FAIL
#SBATCH --account=mh1126
#SBATCH --output=./divfree2d_out.%j.out
#SBATCH --error=./divfree2d_err.%j.out

### ------------------ Input Parameters ---------------- ###
### ------ You MUST edit these lines to set your ------- ###
### ---- build type, directories, the executable(s) ---- ###
### -------- to compile, and your python script -------- ###
### ---------------------------------------------------- ###
buildtype="cuda"
path2CLEO=${HOME}/CLEO/
path2build=${HOME}/CLEO/build_divfree2D/
enableyac=false
executables="divfree2D"

pythonscript=${path2CLEO}/examples/divfreemotion/divfree2d.py
configfile=${path2CLEO}/examples/divfreemotion/src/config/divfree2d_config.yaml
script_args="${configfile}"
### ---------------------------------------------------- ###
### ---------------------------------------------------- ###
### ---------------------------------------------------- ###

### ---------- build, compile and run example ---------- ###
${path2CLEO}/examples/run_example.sh \
  ${buildtype} ${path2CLEO} ${path2build} ${enableyac} \
  "${executables}" ${pythonscript} "${script_args}"
### ---------------------------------------------------- ###
