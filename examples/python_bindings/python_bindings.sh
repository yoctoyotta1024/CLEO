#!/bin/bash
#SBATCH --job-name=pythonbindings
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=10G
#SBATCH --time=00:02:00
#SBATCH --mail-user=clara.bayley@mpimet.mpg.de
#SBATCH --mail-type=FAIL
#SBATCH --account=bm1183
#SBATCH --output=./pythonbindings_out.%j.out
#SBATCH --error=./pythonbindings_err.%j.out

### ---------------------------------------------------- ###
### ------------------ Input Parameters ---------------- ###
### ------ You MUST edit these lines to set your ------- ###
### ---- build type, directories, the executable(s) ---- ###
### -------- to compile, and your python script -------- ###
### ---------------------------------------------------- ###
buildtype="threads"
path2CLEO=${HOME}/CLEO/
path2build=${HOME}/CLEO/build_pybind/
build_flags="-DCLEO_COUPLED_DYNAMICS=numpy -DCLEO_DOMAIN=cartesian \
  -DCLEO_NO_ROUGHPAPER=true -DCLEO_PYTHON=/work/bm1183/m300950/bin/envs/cleoenv/bin/python"
enableyacpython=false
executables="pycleo"

pythonscript=${path2CLEO}/examples/python_bindings/python_bindings.py
configfile=${path2CLEO}/examples/python_bindings/src/config/python_bindings_config.yaml
script_args="${configfile} --do_inputfiles=TRUE --do_run_executable=TRUE --do_plot_results=TRUE"
### ---------------------------------------------------- ###
### ---------------------------------------------------- ###
### ---------------------------------------------------- ###

### ---------- build, compile and run example ---------- ###
${path2CLEO}/examples/run_example_levante.sh \
  ${buildtype} ${path2CLEO} ${path2build} "${build_flags}" ${enableyacpython} \
  "${executables}" ${pythonscript} "${script_args}"
### ---------------------------------------------------- ###
