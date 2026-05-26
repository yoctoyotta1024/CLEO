#!/bin/bash

set -e

### ============================================================ ###
###                      Levante script                          ###
### ============================================================ ###
###
### Usage:
###   ./build_compile_run_plot_cleo.sh [experiment] [buildtype] [compilername] \
###                                    [path2CLEO] [path2build] [build_flags] \
###                                    [yacyaxtroot] [enabledebug] [make_clean] \
###                                    [stacksize_limit]
###
###   All arguments are optional — defaults come from environment variables
###   and example_params.sh.
###
### Arguments:
###   $1  experiment       Name of experiment to run        (default: set below)
###   $2  buildtype        serial | threads | openmp | cuda (default: from experiment)
###   $3  compilername     gcc | intel                      (default: intel)
###   $4  path2CLEO        Absolute path to CLEO source     (default: $CLEO_PATH2CLEO)
###   $5  path2build       Absolute path for build dir      (default: from experiment)
###   $6  build_flags      Extra CMake flags                (default: from experiment)
###   $7  yacyaxtroot      Path to YAC+YAXT install         (default: $CLEO_YACYAXTROOT)
###   $8  enabledebug      true | false                     (default: false)
###   $9  make_clean       true | false                     (default: false)
###   $10 stacksize_limit  ulimit -s value (kB)             (default: 204800)
###
### Available examples (see common/examples/example_params.sh for full details):
###   as2017  cuspbifurc  breakup  shima2009  constthermo2d  divfree2d
###   eurec4a1d  rainshaft1d  python_bindings  kokkostools
###   fromfile  fromfile_irreg  bubble3d
### ============================================================ ###

### ----------------- configure -------------------- ###
export CLEO_MACHINE="levante"
export CLEO_PYTHON=${CLEO_PYTHON:-$(which python3)}

experiment=${1:-"as2017"}
buildtype=${2:-"openmp"}
compilername=${3:-gcc}
path2CLEO=${4:-${HOME}/CLEO}

if [ "${path2CLEO}" == "" ]; then
  echo "Please provide path to CLEO source directory"
  exit 1
fi

[ -f "${path2CLEO}/scripts_2/common/examples/example_params.sh" ] && \
  source ${path2CLEO}/scripts_2/common/examples/example_params.sh "$5" "$6" "${experiment}"

yacyaxtroot=${7:-${HOME}/yacyaxt/${compilername}}
enabledebug=${8:-false}
make_clean=${9:-false}
stacksize_limit=${10:-204800}
### ---------------------------------------------------- ###

### ----------------- export inputs -------------------- ###
export CLEO_BUILDTYPE=${buildtype}
export CLEO_COMPILERNAME=${compilername}
export CLEO_PATH2CLEO=${path2CLEO}
export CLEO_YACYAXTROOT=${yacyaxtroot}
export CLEO_ENABLEDEBUG=${enabledebug}
export CLEO_MAKE_JOBS=32
### ---------------------------------------------------- ###

### -------------------- check inputs ------------------ ###
source ${path2CLEO}/scripts_2/common/bash/src/check_inputs.sh

check_args_not_empty "${CLEO_BUILDTYPE}" "${CLEO_COMPILERNAME}" "${CLEO_PATH2CLEO}" \
                     "${CLEO_PATH2BUILD}" "${CLEO_BUILD_FLAGS}" "${CLEO_YACYAXTROOT}" \
                     "${CLEO_ENABLEDEBUG}" "${stacksize_limit}"
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

### ---------------- compile experiment ---------------- ###
compilecmd="${CLEO_PATH2CLEO}/scripts_2/common/bash/compile_cleo.sh"
echo "${compilecmd} \"${executables}\" ${make_clean}"
source "${compilecmd}" "${executables}" "${make_clean}"
### ---------------------------------------------------- ###

### ------------- load runtime environment ------------- ###
source ${path2CLEO}/scripts_2/levante/bash/src/runtime_settings.sh ${stacksize_limit}
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
