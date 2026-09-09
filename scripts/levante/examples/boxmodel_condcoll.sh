#!/bin/bash
#SBATCH --job-name=boxmodel_condcoll
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --gpus=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=10G
#SBATCH --time=00:20:00
#SBATCH --mail-user=clara.bayley@mpimet.mpg.de
#SBATCH --mail-type=FAIL
#SBATCH --account=mh0731
#SBATCH --output=./boxmodel_out.%j.out
#SBATCH --error=./boxmodel_err.%j.out

### ---------------------------------------------------- ###
### ------------------ Input Parameters ---------------- ###
### ------ You MUST edit these lines to set your ------- ###
### ---- build type, directories, the executable(s) ---- ###
### -------- to compile, and your python script -------- ###
### ---------------------------------------------------- ###
do_build="true"
buildtype="openmp"
compilername="intel"
path2CLEO=${HOME}/CLEO/
path2build=${HOME}/CLEO/build_colls0d/boxmodel_condcoll/
build_flags="-DCLEO_COUPLED_DYNAMICS=fromfile -DCLEO_DOMAIN=cartesian \
  -DCLEO_NO_ROUGHPAPER=true -DCLEO_NO_PYBINDINGS=true"
executables="boxmodel_condcoll"

pythonscript=${path2CLEO}/examples/boxmodel_condcoll/boxmodel_condcoll.py
src_config_filename=${path2CLEO}/examples/boxmodel_condcoll/src/config/boxmodel_condcoll_config.yaml
script_args="${src_config_filename} --do_inputfiles --do_run_executable --do_plot_results"
### ---------------------------------------------------- ###
### ---------------------------------------------------- ###
### ---------------------------------------------------- ###

### ---------- build, compile and run example ---------- ###
${path2CLEO}/scripts/levante/examples/build_compile_run_plot.sh ${do_build} \
  ${buildtype} ${compilername} ${path2CLEO} ${path2build} "${build_flags}" \
  "${executables}" ${pythonscript} "${script_args}"
### ---------------------------------------------------- ###
