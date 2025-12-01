#!/bin/bash

### ------ Generic script to build CLEO, compile some ------ ###
### ----- of its executables and run a python scripts  ----- ###

### ---------------------------------------------------- ###
### ------------------ Input Parameters ---------------- ###
### ------ You MUST edit these lines to set your ------- ###
### ----- environment, build type, directories, the ---- ###
### --------- executable(s) to compile and your -------- ###
### --------------  python script to run. -------------- ###
### ---------------------------------------------------- ###
do_build=$1  # == "true" or otherwise false
buildtype=$2
compilername=$3
path2CLEO=$4
path2build=$5
build_flags=$6
executables="$7"
pythonscript=$8
script_args="$9"

python=${CLEO_PYTHON}
yacyaxtroot=${CLEO_YACYAXTROOT}

enabledebug=false
make_clean=false
### ---------------------------------------------------- ###
### ---------------------------------------------------- ###
### ---------------------------------------------------- ###

### -------------------- print inputs ------------------ ###
echo "----- Running Example -----"
echo "buildtype = ${buildtype}"
echo "compilername = ${compilername}"
echo "path2CLEO = ${path2CLEO}"
echo "path2build = ${path2build}"
echo "build_flags = ${build_flags}"
echo "yacyaxtroot = ${yacyaxtroot}"
echo "executables = ${executables}"
echo "pythonscript = ${pythonscript}"
echo "script_args = ${script_args}"
echo "---------------------------"
### ---------------------------------------------------- ###

### --------------- build and compile CLEO ------------- ###
if [ "${do_build}" == "true" ]
then
  cmd="${path2CLEO}/scripts/vanilla/build_compile_cleo.sh \
    ${buildtype}
    ${compilername}
    ${path2CLEO}
    ${path2build}
    ${yacyaxtroot}
    "\"${build_flags}\""
    "\"${executables}\""
    ${enabledebug}
    ${make_clean}"
  echo ${cmd}
  eval ${cmd}
fi
### ---------------------------------------------------- ###

### --------- run model through Python script ---------- ###
export CLEO_PATH2CLEO=${path2CLEO}
export CLEO_BUILDTYPE=${buildtype}
export CLEO_COMPILERNAME=${compilername}
export CLEO_YACYAXTROOT=${yacyaxtroot}
source ${path2CLEO}/scripts/vanilla/bash/src/runtime_settings.sh

# TODO(ALL): split python scripts away from running executable
${python} ${pythonscript} ${path2CLEO} ${path2build} ${script_args}
### ---------------------------------------------------- ###
