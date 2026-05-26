#!/bin/bash

set -e

### ============================================================ ###
###                      Vanilla script                          ###
### ============================================================ ###
###
### Usage:
###   ./build_run_cleo.sh [experiment] [buildtype] [compilername] \
###                       [path2CLEO] [path2build] [yacyaxtroot] [build_flags] \
###                       [enabledebug] [make_clean]
###
###   All arguments are optional — defaults come from environment variables
###   and example_params.sh.
###
### Arguments:
###   $1  experiment     Name of experiment to run       (default: set below)
###   $2  buildtype      serial | threads | openmp       (default: from experiment)
###   $3  compilername   gcc                             (default: from experiment)
###   $4  path2CLEO      Absolute path to CLEO source    (default: $CLEO_PATH2CLEO)
###   $5  path2build     Absolute path for build dir     (default: build_<experiment>)
###   $6  build_flags    Extra CMake flags                (default: from experiment)
###   $7  yacyaxtroot    Path to YAC+YAXT install        (default: $CLEO_YACYAXTROOT)
###   $8  enabledebug    true | false                    (default: false)
###   $9  make_clean     true | false                    (default: true)
###
### Available examples (see common/examples/example_params.sh for full details):
###   as2017  cuspbifurc  breakup  shima2009  constthermo2d  divfree2d
###   eurec4a1d  rainshaft1d  python_bindings  kokkostools
###   fromfile  fromfile_irreg  bubble3d
### ============================================================ ###

#step 1: configure
export CLEO_MACHINE="vanilla"

experiment=${1:-"as2017"}
buildtype=${2:-"serial"}
compilername=${3:-gcc}
path2CLEO=${4:-${CLEO_PATH2CLEO}}

if [ "${path2CLEO}" == "" ]; then
  echo "Please provide path to CLEO source directory"
  exit 1
fi

[ -f "${path2CLEO}/scripts_2/common/examples/example_params.sh" ] && \
  source ${path2CLEO}/scripts_2/common/examples/example_params.sh  "$5" "$6" "${experiment}"

yacyaxtroot=${7:-${CLEO_YACYAXTROOT}}
enabledebug=${8:-false}
make_clean=${9:-false}

### ----------------- export inputs -------------------- ###
export CLEO_BUILDTYPE=${buildtype}
export CLEO_COMPILERNAME=${compilername}
export CLEO_PATH2CLEO=${path2CLEO}
export CLEO_YACYAXTROOT=${yacyaxtroot}
export CLEO_ENABLEDEBUG=${enabledebug}
### ---------------------------------------------------- ###

### -------------------- check inputs ------------------ ###
source ${path2CLEO}/scripts_2/common/bash/src/check_inputs.sh

check_args_not_empty "${CLEO_BUILDTYPE}" "${CLEO_COMPILERNAME}" "${CLEO_PATH2CLEO}" \
                     "${CLEO_PATH2BUILD}" "${CLEO_BUILD_FLAGS}" "${CLEO_YACYAXTROOT}" \
                     "${CLEO_ENABLEDEBUG}"

### ---------------------------------------------------- ###

[ -f "${path2CLEO}/scripts_2/common/bash/src/print_configuration.sh" ] && \
  source "${path2CLEO}/scripts_2/common/bash/src/print_configuration.sh" "${experiment}"

### --------------------- build CLEO ------------------- ###
buildcmd="${CLEO_PATH2CLEO}/scripts_2/common/bash/build_cleo.sh"
if [ ! -f "${buildcmd}" ]; then
  echo "Error: build script not found at ${buildcmd}"
  exit 1
fi
echo "${buildcmd}"
source "${buildcmd}"
### ---------------------------------------------------- ###

### ---------------- compile experiment --------------- ###
compilecmd="${CLEO_PATH2CLEO}/scripts_2/common/bash/compile_cleo.sh"
echo "${compilecmd} \"${executables}\" ${make_clean}"
source "${compilecmd}" "${executables}" "${make_clean}"
### ---------------------------------------------------- ###

### ----------- run Python plot/analysis script -------- ###
if [ -z "${pythonscript}" ]; then
  echo "Error: no Python script provided."
  exit 1
fi
if [ ! -f "${pythonscript}" ]; then
  echo "Error: Python script not found at ${pythonscript}"
  exit 1
fi
if [ -z "${CLEO_PYTHON}" ]; then
  echo "Error: CLEO_PYTHON is not set."
  exit 1
fi

echo "Running: ${CLEO_PYTHON} ${pythonscript} ${path2CLEO} ${CLEO_PATH2BUILD} ${script_args}"
${CLEO_PYTHON} "${pythonscript}" "${path2CLEO}" "${CLEO_PATH2BUILD}" ${script_args}
### ---------------------------------------------------- ###
