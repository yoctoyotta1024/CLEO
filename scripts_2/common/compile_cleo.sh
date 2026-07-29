#!/bin/bash

### Core make logic for compiling CLEO executables.
### Requires CLEO_* environment variables to already be exported.
### Machine-specific module loading must be done by the calling wrapper BEFORE
### sourcing this script.

set -e

executables=$1   # space-separated list, or "NONE" to build all libraries
make_clean=$2
make_jobs=${CLEO_MAKE_JOBS:-8}  # override by setting CLEO_MAKE_JOBS (e.g. 128 on Levante)

### ---------------------------------------------------- ###

### ---------------- compile executables --------------- ###

cd ${CLEO_PATH2BUILD} && pwd

if [ "${make_clean}" == "true" ]; then
  cmd="make clean"
  echo ${cmd}
  eval ${cmd}
fi

if [ "${executables}" == "NONE" ]; then
  cmd="make -j ${make_jobs}"
else
  cmd="make -j ${make_jobs} ${executables}"
fi
echo ${cmd}
eval ${cmd}
### ---------------------------------------------------- ###
