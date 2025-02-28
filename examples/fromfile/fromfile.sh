#!/bin/bash
#SBATCH --job-name=fromfile
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=16
#SBATCH --mem=10G
#SBATCH --time=00:05:00
#SBATCH --mail-user=clara.bayley@mpimet.mpg.de
#SBATCH --mail-type=FAIL
#SBATCH --account=bm1183
#SBATCH --output=./fromfile_out.%j.out
#SBATCH --error=./fromfile_err.%j.out

### ---------------------------------------------------- ###
### ------------------ Input Parameters ---------------- ###
### ------ You MUST edit these lines to set your ------- ###
### ---- build type, directories, the executable(s) ---- ###
### -------- to compile, and your python script -------- ###
### ---------------------------------------------------- ###
buildtype="openmp"
path2CLEO=${HOME}/CLEO/
path2build=${HOME}/CLEO/build_fromfile/
build_flags="-DCLEO_COUPLED_DYNAMICS=fromfile -DCLEO_DOMAIN=cartesian -DCLEO_NO_ROUGHPAPER=true"
enableyac=false
executables="fromfile"

pythonscript=${path2CLEO}/examples/fromfile/fromfile.py
configfile=${path2CLEO}/examples/fromfile/src/config/fromfile_config.yaml
script_args="${configfile} --do_inputfiles=TRUE --do_run_executable=TRUE --do_plot_results=TRUE --ntasks=4"
### ---------------------------------------------------- ###
### ---------------------------------------------------- ###
### ---------------------------------------------------- ###

### ---------- build, compile and run example ---------- ###
${path2CLEO}/examples/run_example_levante.sh \
  ${buildtype} ${path2CLEO} ${path2build} "${build_flags}" ${enableyac} \
  "${executables}" ${pythonscript} "${script_args}"
### ---------------------------------------------------- ###
