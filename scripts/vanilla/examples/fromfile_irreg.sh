#!/bin/bash

### ---------------------------------------------------- ###
### ------------------ Input Parameters ---------------- ###
### ------ You MUST edit these lines to set your ------- ###
### ---- build type, directories, the executable(s) ---- ###
### -------- to compile, and your python script -------- ###
### ---------------------------------------------------- ###
do_build="true"
buildtype="openmp"
compilername="gcc"
path2CLEO=${CLEO_PATH2CLEO}
path2build=${path2CLEO}/build_fromfile_irreg/
build_flags="-DCLEO_COUPLED_DYNAMICS=fromfile -DCLEO_DOMAIN=cartesian \
  -DCLEO_NO_ROUGHPAPER=true -DCLEO_NO_PYBINDINGS=true"
executables="fromfile_irreg"

pythonscript=${path2CLEO}/examples/fromfile_irreg/fromfile_irreg.py
src_config_filename=${path2CLEO}/examples/fromfile_irreg/src/config/fromfile_irreg_config.yaml
script_args="${src_config_filename} \
  --do_inputfiles --do_run_executable --do_plot_results --ntasks=4"
### ---------------------------------------------------- ###
### ---------------------------------------------------- ###
### ---------------------------------------------------- ###

### ---------- build, compile and run example ---------- ###
${path2CLEO}/scripts/vanilla/examples/build_compile_run_plot.sh ${do_build} \
  ${buildtype} ${compilername} ${path2CLEO} ${path2build} "${build_flags}" \
  "${executables}" ${pythonscript} "${script_args}"
### ---------------------------------------------------- ###
