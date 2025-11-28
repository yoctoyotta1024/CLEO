#!/bin/bash

### Please note: script may assume required CLEO_[XXX]
### variables have already exported (!)

set -e

executables=$1                     # if == "NONE" only libraries built
make_clean=$2

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
bashsrc=${SCRIPT_DIR}/src

### -------------------- check inputs ------------------ ###
source ${bashsrc}/check_inputs.sh
check_args_not_empty "${CLEO_BUILDTYPE}" "${CLEO_COMPILERNAME}" \
  "${CLEO_PATH2CLEO}" "${CLEO_PATH2BUILD}"

check_source_and_build_paths
check_buildtype
check_compilername
### ---------------------------------------------------- ###

### ---------------- compile executables --------------- ###
echo "### --------------- Compile Inputs -------------- ###"
echo "CLEO_BUILDTYPE: ${CLEO_BUILDTYPE}"
echo "CLEO_COMPILERNAME: ${CLEO_COMPILERNAME}"
echo "CLEO_PATH2CLEO: ${CLEO_PATH2CLEO}"
echo "CLEO_PATH2BUILD: ${CLEO_PATH2BUILD}"

echo "executables: ${executables}"
echo "make_clean: ${make_clean}"
echo "### ------------------------------------------- ###"

cd ${CLEO_PATH2BUILD} && pwd

if [ "${make_clean}" == "true" ]
then
  cmd="make clean"
  echo ${cmd}
  eval ${cmd}
fi

if [ ${executables} == "NONE" ]
then
  cmd="make -j 128"
else
  cmd="make -j 128 ${executables}"
fi
echo ${cmd}
eval ${cmd}
### ---------------------------------------------------- ###
