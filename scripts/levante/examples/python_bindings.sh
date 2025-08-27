#!/bin/bash
#SBATCH --job-name=pythonbindings
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=10G
#SBATCH --time=00:05:00
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
do_build="true"
buildtype="threads"
compilername="gcc"
path2CLEO=${HOME}/CLEO/
path2build=${HOME}/CLEO/build_pybind/
build_flags="-DCLEO_COUPLED_DYNAMICS=numpy -DCLEO_DOMAIN=cartesian \
  -DCLEO_NO_ROUGHPAPER=true -DCLEO_PYTHON=/home/m/m300950/CLEO/.venv/bin/python3"
executables="cleo_python_bindings"

pythonscript=${path2CLEO}/examples/python_bindings/python_bindings.py
src_config_filename=${path2CLEO}/examples/python_bindings/src/config/pybind_config.yaml
script_args="${src_config_filename} \
  --do_inputfiles --do_run_executable --do_plot_results"
### ---------------------------------------------------- ###
### ---------------------------------------------------- ###
### ---------------------------------------------------- ###

### ---------- build, compile and run example ---------- ###
${path2CLEO}/scripts/levante/examples/build_compile_run_plot.sh ${do_build} \
  ${buildtype} ${compilername} ${path2CLEO} ${path2build} "${build_flags}" \
  "${executables}" ${pythonscript} "${script_args}"
### ---------------------------------------------------- ###
