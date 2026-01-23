#!/bin/bash

set -e

### ------------------ input parameters ---------------- ###
### ----- You need to edit these lines to specify ------ ###
### ----- your build configuration and executables ----- ###
### ---------------------------------------------------- ###
buildtype=$1                                   # "serial", "threads", or "openmp"
compilername=${2:-gcc}                         # "gcc"
path2CLEO=${3:-${CLEO_PATH2CLEO}}              # must be absolute path
path2build=${4:-${path2CLEO}/build}            # should be absolute path
yacyaxtroot=${5:-${CLEO_YACYAXTROOT}}          # yac and yaxt in yacyaxtroot/yac and yacyaxtroot/yaxt
build_flags=${6:-"-DCLEO_COUPLED_DYNAMICS="" -DCLEO_PYTHON=${CLEO_PYTHON}"} # CLEO_BUILD_FLAGS
executables=${7:-"cleocoupledsdm"}             # list of executables to compile or "NONE"
enabledebug=${8:-false}                        # == "true" or otherwise false
make_clean=${9:-true}                          # == "true" or otherwise false
### ---------------------------------------------------- ###

### ------------------ check arguments ----------------- ###
if [ "${path2CLEO}" == "" ]
then
  echo "Please provide path to CLEO source directory"
  exit 1
fi
source ${path2CLEO}/scripts/vanilla/bash/src/check_inputs.sh
check_args_not_empty "${buildtype}" "${path2CLEO}" "${path2build}"
check_args_not_empty "${compilername}" "${yacyaxtroot}" "${enabledebug}"
### ---------------------------------------------------- ###

### ----------------- export inputs -------------------- ###
export CLEO_BUILDTYPE=${buildtype}
export CLEO_COMPILERNAME=${compilername}
export CLEO_PATH2CLEO=${path2CLEO}
export CLEO_PATH2BUILD=${path2build}
export CLEO_BUILD_FLAGS=${build_flags}
export CLEO_YACYAXTROOT=${yacyaxtroot}
export CLEO_ENABLEDEBUG=${enabledebug}
### ---------------------------------------------------- ###

### -------------------- check inputs ------------------ ###
check_source_and_build_paths
check_buildtype
check_compilername
check_yac
### ---------------------------------------------------- ###

### -------------------- print inputs ------------------- ###
echo "### --------------- User Inputs -------------- ###"
echo "CLEO_BUILDTYPE = ${CLEO_BUILDTYPE}"
echo "CLEO_COMPILERNAME = ${CLEO_COMPILERNAME}"
echo "CLEO_PATH2CLEO = ${CLEO_PATH2CLEO}"
echo "CLEO_PATH2BUILD = ${CLEO_PATH2BUILD}"
echo "CLEO_BUILD_FLAGS = ${CLEO_BUILD_FLAGS}"
echo "CLEO_YACYAXTROOT = ${CLEO_YACYAXTROOT}"
echo "CLEO_ENABLEDEBUG = ${CLEO_ENABLEDEBUG}"
echo "executables = ${executables}"
echo "### ------------------------------------------- ###"
### ---------------------------------------------------- ###

### --------------------- build CLEO ------------------- ###
buildcmd="${CLEO_PATH2CLEO}/scripts/vanilla/bash/build_cleo.sh"
echo ${buildcmd}
eval ${buildcmd}
### ---------------------------------------------------- ###

### ---------------- compile executables --------------- ###
compilecmd="${CLEO_PATH2CLEO}/scripts/vanilla/bash/compile_cleo.sh \"${executables}\" ${make_clean}"
echo ${compilecmd}
eval ${compilecmd}
### ---------------------------------------------------- ###
