#!/bin/bash
#SBATCH --job-name=shima2009
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --gpus=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=10G
#SBATCH --time=00:20:00
#SBATCH --mail-user=clara.bayley@mpimet.mpg.de
#SBATCH --mail-type=FAIL
#SBATCH --account=bm1183
#SBATCH --output=./shima2009_out.%j.out
#SBATCH --error=./shima2009_err.%j.out

### ---------------------------------------------------- ###
### ------------------ Input Parameters ---------------- ###
### ------ You MUST edit these lines to set your ------- ###
### ---- build type, directories, the executable(s) ---- ###
### -------- to compile, and your python script -------- ###
### ---------------------------------------------------- ###
buildtype="cuda"
path2CLEO=${HOME}/CLEO/
path2build=${HOME}/CLEO/build_colls0d/${buildtype}/
build_flags="-DCLEO_COUPLED_DYNAMICS=null -DCLEO_DOMAIN=cartesian -DCLEO_NO_ROUGHPAPER=true"
enableyac=false
executables="golcolls longcolls"

pythonscript=${path2CLEO}/examples/boxmodelcollisions/shima2009.py
configfile=${path2CLEO}/examples/boxmodelcollisions/shima2009_config.yaml
script_args="${configfile} golovin long1 long2"
### ---------------------------------------------------- ###
### ---------------------------------------------------- ###
### ---------------------------------------------------- ###

### ---------- build, compile and run example ---------- ###
${path2CLEO}/examples/run_example_levante.sh \
  ${buildtype} ${path2CLEO} ${path2build} "${build_flags}" ${enableyac} \
  "${executables}" ${pythonscript} "${script_args}"
### ---------------------------------------------------- ###
