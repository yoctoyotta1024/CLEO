#!/bin/bash
#SBATCH --job-name=bubble3d
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=256
#SBATCH --mem=10G
#SBATCH --time=00:05:00
#SBATCH --mail-user=clara.bayley@mpimet.mpg.de
#SBATCH --mail-type=FAIL
#SBATCH --account=bm1183
#SBATCH --output=./bubble3d_out.%j.out
#SBATCH --error=./bubble3d_err.%j.out

# TODO(all): complete python script(s) for example (and fix MPI linker error?)

### ---------------------------------------------------- ###
### ------------------ Input Parameters ---------------- ###
### ------ You MUST edit these lines to set your ------- ###
### ---- build type, directories, the executable(s) ---- ###
### -------- to compile, and your python script -------- ###
### ---------------------------------------------------- ###
buildtype="openmp"
path2CLEO=${HOME}/CLEO/
path2build=${HOME}/CLEO/build_bubble3d/
build_flags="-DCLEO_COUPLED_DYNAMICS=yac -DCLEO_DOMAIN=cartesian \
  -DCLEO_NO_ROUGHPAPER=true -DCLEO_NO_PYBINDINGS=true"
executables="bubble3d"

pythonscript=${path2CLEO}/examples/bubble3d/bubble3d.py
configfile=${path2CLEO}/examples/bubble3d/src/config/bubble3d_config.yaml
script_args="${configfile}"
### ---------------------------------------------------- ###
### ---------------------------------------------------- ###
### ---------------------------------------------------- ###

### ---------- build, compile and run example ---------- ###
${path2CLEO}/examples/run_example_levante.sh \
  ${buildtype} ${path2CLEO} ${path2build} "${build_flags}" \
  "${executables}" ${pythonscript} "${script_args}"
### ---------------------------------------------------- ###
