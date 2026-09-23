#!/bin/bash

### Usage: print_configuration <experiment>
###
### Prints the current build configuration to stdout.

print_configuration() {
  local experiment_name=$1

  echo "### --------------- User Inputs -------------- ###"
  echo "CLEO_MACHINE = ${CLEO_MACHINE}"
  echo "EXPERIMENT = ${experiment_name}"
  echo "CLEO_BUILDTYPE = ${CLEO_BUILDTYPE}"
  echo "CLEO_COMPILERNAME = ${CLEO_COMPILERNAME}"
  echo "CLEO_PATH2CLEO = ${CLEO_PATH2CLEO}"
  echo "CLEO_PATH2BUILD = ${CLEO_PATH2BUILD}"
  echo "CLEO_BUILD_FLAGS = ${CLEO_BUILD_FLAGS}"
  echo "CLEO_YACYAXTROOT = ${CLEO_YACYAXTROOT}"
  echo "CLEO_ENABLEDEBUG = ${CLEO_ENABLEDEBUG}"
  echo "CLEO_PYTHON = ${CLEO_PYTHON}"
  echo "CLEO_MAKE_JOBS = ${CLEO_MAKE_JOBS}"
  echo "STACKSIZE_LIMIT = ${stacksize_limit:-(unchanged)}"
  echo "STEPS = ${steps:-all}"
  echo "### ------------------------------------------- ###"
}

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  print_configuration "$@"
fi
