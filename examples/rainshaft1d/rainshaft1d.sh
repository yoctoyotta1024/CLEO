#!/bin/bash
#SBATCH --job-name=rain1d
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=256
#SBATCH --mem=10G
#SBATCH --time=00:10:00
#SBATCH --mail-user=clara.bayley@mpimet.mpg.de
#SBATCH --mail-type=FAIL
#SBATCH --account=bm1183
#SBATCH --output=./rain1d_out.%j.out
#SBATCH --error=./rain1d_err.%j.out

### ------------------ Input Parameters ---------------- ###
### ------ You MUST edit these lines to set your ------- ###
### ---- build type, directories, the executable(s) ---- ###
### -------- to compile, and your python script -------- ###
### ---------------------------------------------------- ###
buildtype="openmp"
path2CLEO=${HOME}/CLEO/
path2build=${HOME}/CLEO/build_rshaft1d/
build_flags="-DCLEO_COUPLED_DYNAMICS=fromfile -DCLEO_DOMAIN=cartesian \
  -DCLEO_NO_ROUGHPAPER=true -DCLEO_NO_PYBINDINGS=true"
executables="rshaft1d"

pythonscript=${path2CLEO}/examples/rainshaft1d/rainshaft1d.py
configfile=${path2CLEO}/examples/rainshaft1d/src/config/rain1d_config.yaml
script_args="${configfile}"
### ---------------------------------------------------- ###
### ---------------------------------------------------- ###
### ---------------------------------------------------- ###

### ---------- build, compile and run example ---------- ###
${path2CLEO}/examples/run_example_levante.sh \
  ${buildtype} ${path2CLEO} ${path2build} "${build_flags}" \
  "${executables}" ${pythonscript} "${script_args}"
### ---------------------------------------------------- ###
