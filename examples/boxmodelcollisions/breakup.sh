#!/bin/bash
#SBATCH --job-name=breakup
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --gpus=4
#SBATCH --ntasks-per-node=128
#SBATCH --mem=30G
#SBATCH --time=00:10:00
#SBATCH --mail-user=clara.bayley@mpimet.mpg.de
#SBATCH --mail-type=FAIL
#SBATCH --account=mh1126
#SBATCH --output=./breakup_out.%j.out
#SBATCH --error=./breakup_err.%j.out

### ---------------------------------------------------- ###
### ------------------ Input Parameters ---------------- ###
### ------ You MUST edit these lines to set your ------- ###
### ---- build type, directories, the executable(s) ---- ###
### -------- to compile, and your python script -------- ###
### ---------------------------------------------------- ###
buildtype="cuda"
path2CLEO=${HOME}/CLEO/
path2build=${HOME}/CLEO/build_colls0D/
enableyac=false
executables="longcolls lowlistcolls szakallurbichcolls testikstraubcolls"

pythonscript=${path2CLEO}/examples/boxmodelcollisions/breakup.py
configfile=${path2CLEO}/examples/boxmodelcollisions/breakup_config.yaml
script_args="${configfile} long lowlist szakallurbich testikstraub"
### ---------------------------------------------------- ###
### ---------------------------------------------------- ###
### ---------------------------------------------------- ###

### ---------- build, compile and run example ---------- ###
${path2CLEO}/examples/run_example.sh \
  ${buildtype} ${path2CLEO} ${path2build} ${enableyac} \
  "${executables}" ${pythonscript} "${script_args}"
### ---------------------------------------------------- ###
