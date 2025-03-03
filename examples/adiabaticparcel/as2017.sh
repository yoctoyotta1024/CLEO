#!/bin/bash
#SBATCH --job-name=as2017
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=10G
#SBATCH --time=00:10:00
#SBATCH --mail-user=clara.bayley@mpimet.mpg.de
#SBATCH --mail-type=FAIL
#SBATCH --account=bm1183
#SBATCH --output=./as2017_out.%j.out
#SBATCH --error=./as2017_err.%j.out

### ---------------------------------------------------- ###
### ------------------ Input Parameters ---------------- ###
### ------ You MUST edit these lines to set your ------- ###
### ---- build type, directories, the executable(s) ---- ###
### -------- to compile, and your python script -------- ###
### ---------------------------------------------------- ###
buildtype="openmp"
path2CLEO=${HOME}/CLEO/
path2build=${HOME}/CLEO/build_adia0d/${buildtype}/
build_flags="-DCLEO_COUPLED_DYNAMICS=cvode -DCLEO_DOMAIN=cartesian -DCLEO_NO_ROUGHPAPER=true"
enableyac=false
executables="adia0d"

pythonscript=${path2CLEO}/examples/adiabaticparcel/as2017.py
configfile=${path2CLEO}/examples/adiabaticparcel/src/config/as2017_config.yaml
script_args="${configfile}"
### ---------------------------------------------------- ###
### ---------------------------------------------------- ###
### ---------------------------------------------------- ###

### ---------- build, compile and run example ---------- ###
${path2CLEO}/examples/run_example_levante.sh \
  ${buildtype} ${path2CLEO} ${path2build} "${build_flags}" ${enableyac} \
  "${executables}" ${pythonscript} "${script_args}"
### ---------------------------------------------------- ###
