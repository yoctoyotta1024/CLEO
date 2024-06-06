#!/bin/bash
#SBATCH --job-name=yac_3d
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --mem=30G
#SBATCH --time=00:05:00
#SBATCH --mail-user=clara.bayley@mpimet.mpg.de
#SBATCH --mail-type=FAIL
#SBATCH --account=mh1126
#SBATCH --output=./yac_3d_out.%j.out
#SBATCH --error=./yac_3d_err.%j.out

# TODO(all): python script(s) for example and fix MPI linker error

### ---------------------------------------------------- ###
### ------------------ Input Parameters ---------------- ###
### ------ You MUST edit these lines to set your ------- ###
### ---- build type, directories, the executable(s) ---- ###
### -------- to compile, and your python script -------- ###
### ---------------------------------------------------- ###
buildtype="openmp"
path2CLEO=${HOME}/CLEO/
path2build=${HOME}/CLEO/build_yac3d/
enableyac=true
executables="yac_3d"

pythonscript=${path2CLEO}/examples/yac/yac_3d/yac_3d_fromfile.py
configfile=${path2CLEO}/examples/yac/yac_3d/src/config/yac1_fromfile_config.yaml
script_args="${configfile}"
### ---------------------------------------------------- ###
### ---------------------------------------------------- ###
### ---------------------------------------------------- ###

### ---------- build, compile and run example ---------- ###
${path2CLEO}/examples/run_example.sh \
  ${buildtype} ${path2CLEO} ${path2build} ${enableyac} \
  "${executables}" ${pythonscript} "${script_args}"
### ---------------------------------------------------- ###
