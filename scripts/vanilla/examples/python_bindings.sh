#!/bin/bash

### ---------------------------------------------------- ###
### ------------------ Input Parameters ---------------- ###
### ------ You MUST edit these lines to set your ------- ###
### ---- build type, directories, the executable(s) ---- ###
### -------- to compile, and your python script -------- ###
### ---------------------------------------------------- ###
do_build="true"
buildtype="threads"
compilername="gcc"
path2CLEO=${CLEO_PATH2CLEO}
path2build=${path2CLEO}/build_pybind/
build_flags="-DCLEO_COUPLED_DYNAMICS=numpy -DCLEO_DOMAIN=cartesian \
  -DCLEO_NO_ROUGHPAPER=true -DCLEO_PYTHON=${CLEO_PYTHON}"
executables="cleo_python_bindings"

pythonscript=${path2CLEO}/examples/python_bindings/python_bindings.py
src_config_filename=${path2CLEO}/examples/python_bindings/src/config/pybind_config.yaml
script_args="${src_config_filename} \
  --do_inputfiles --do_run_executable --do_plot_results"
### ---------------------------------------------------- ###
### ---------------------------------------------------- ###
### ---------------------------------------------------- ###

if [[ "${compilername}" != "gcc" ]]
then
  echo "python bindings example currently only working with gcc compiler"
  echo "-> please use compilername=gcc"
  exit 1
fi

### ---------- build, compile and run example ---------- ###
${path2CLEO}/scripts/vanilla/examples/build_compile_run_plot.sh ${do_build} \
  ${buildtype} ${compilername} ${path2CLEO} ${path2build} "${build_flags}" \
  "${executables}" ${pythonscript} "${script_args}"
### ---------------------------------------------------- ###
