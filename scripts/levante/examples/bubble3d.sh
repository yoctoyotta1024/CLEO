#!/bin/bash
#SBATCH --job-name=bubble3d
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=32
#SBATCH --mem=10G
#SBATCH --time=00:05:00
#SBATCH --mail-user=clara.bayley@mpimet.mpg.de
#SBATCH --mail-type=FAIL
#SBATCH --account=mh0731
#SBATCH --output=./bubble3d_out.%j.out
#SBATCH --error=./bubble3d_err.%j.out

### TODO(CB): vanilla script for bubble example

### ---------------------------------------------------- ###
### ------------------ Input Parameters ---------------- ###
### ------ You MUST edit these lines to set your ------- ###
### ---- build type, directories, the executable(s) ---- ###
### -------- to compile, and your python script -------- ###
### ---------------------------------------------------- ###
do_build="true"
buildtype="openmp"
compilername="gcc"
path2CLEO=${HOME}/CLEO/
path2build=/work/mh0731/m300950/icon-mpim/build_cleo
# path2build=${HOME}/CLEO/build_bubble3d/
path2experiment="/work/mh0731/m300950/icon-mpim/experiments/bubble_cleo"
path2iconfiles="/work/mh0731/m300950/icon-mpim/experiments/bubble_1mom/outdata"
build_flags="-DCLEO_COUPLED_DYNAMICS=yac -DCLEO_DOMAIN=cartesian \
  -DCLEO_NO_ROUGHPAPER=true -DCLEO_NO_PYBINDINGS=true"
executables="bubble3d"

pythonscript=${path2CLEO}/examples/bubble3d/bubble3d.py
src_config_filename=${path2CLEO}/examples/bubble3d/src/config/bubble3d_config.yaml
src_yac_config_filename=${path2CLEO}/examples/bubble3d/src/config/yac_icon_cleo_coupling_config.yaml
script_args="${src_config_filename} ${src_yac_config_filename} ${path2experiment} ${path2iconfiles} \
  --do_inputfiles" # --do_plot_results"
### ---------------------------------------------------- ###
### ---------------------------------------------------- ###
### ---------------------------------------------------- ###

if [[ "${compilername}" != "gcc" ]]
then
  echo "bubble3d example currently only working on Levante with gcc compiler"
  echo "-> please use compilername=gcc"
  exit 1
fi

### ---------- build, compile and run example ---------- ###
${path2CLEO}/scripts/levante/examples/build_compile_run_plot.sh ${do_build} \
  ${buildtype} ${compilername} ${path2CLEO} ${path2build} "${build_flags}" \
  "${executables}" ${pythonscript} "${script_args}"
### ---------------------------------------------------- ###
