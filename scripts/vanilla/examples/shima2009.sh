#!/bin/bash

### ---------------------------------------------------- ###
### ------------------ Input Parameters ---------------- ###
### ------ You MUST edit these lines to set your ------- ###
### ---- build type, directories, the executable(s) ---- ###
### -------- to compile, and your python script -------- ###
### ---------------------------------------------------- ###
do_build="true"
buildtype="cuda"
compilername="gcc"
path2CLEO=${HOME}/CLEO/
path2build=${HOME}/CLEO/build_colls0d/shima2009/
build_flags="-DCLEO_COUPLED_DYNAMICS=null -DCLEO_DOMAIN=cartesian \
  -DCLEO_NO_ROUGHPAPER=true -DCLEO_NO_PYBINDINGS=true"
executables="golcolls longcolls"

pythonscript=${path2CLEO}/examples/boxmodelcollisions/shima2009.py
src_config_filename=${path2CLEO}/examples/boxmodelcollisions/src/config/shima2009_config.yaml
script_args="${src_config_filename} --kernels golovin long1 long2 \
  --do_inputfiles --do_run_executable --do_plot_results"
### ---------------------------------------------------- ###
### ---------------------------------------------------- ###
### ---------------------------------------------------- ###

### ---------- build, compile and run example ---------- ###
${path2CLEO}/scripts/levante/examples/build_compile_run_plot.sh ${do_build} \
  ${buildtype} ${compilername} ${path2CLEO} ${path2build} "${build_flags}" \
  "${executables}" ${pythonscript} "${script_args}"
### ---------------------------------------------------- ###
