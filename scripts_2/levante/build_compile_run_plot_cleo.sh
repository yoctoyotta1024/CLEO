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
###   and experiments.sh.
###
### Arguments:
###   $1  experiment       Name of experiment to run        (default: set below)
###   $2  buildtype        serial | threads | openmp | cuda (default: openmp)
###   $3  compilername     gcc | intel                      (default: gcc)
###   $4  path2CLEO        Absolute path to CLEO source     (default: $CLEO_PATH2CLEO)
###   $5  path2build       Absolute path for build dir      (default: from experiment)
###   $6  build_flags      Extra CMake flags                (default: from experiment)
###   $7  yacyaxtroot      Path to YAC+YAXT install         (default: $CLEO_YACYAXTROOT)
###   $8  enabledebug      true | false                     (default: false)
###   $9  make_clean       true | false                     (default: false)
###   $10 stacksize_limit  ulimit -s value (kB)             (default: 204800)
###
### Available examples (see common/experiments.sh for full details):
###   as2017  cuspbifurc  breakup  shima2009  constthermo2d  divfree2d
###   eurec4a1d  rainshaft1d  python_bindings
###   fromfile  fromfile_irreg  bubble3d
### ============================================================ ###

### ----------------- configure -------------------- ###
export CLEO_MACHINE="levante"
export CLEO_PYTHON=${CLEO_PYTHON:-$(which python3)}

experiment=${1:-"as2017"}
buildtype=${2:-"openmp"}
compilername=${3:-gcc}
path2CLEO=${4:-${HOME}/CLEO}

# Export build context early so experiment resolution can use it.
export CLEO_BUILDTYPE=${buildtype}
export CLEO_COMPILERNAME=${compilername}
export CLEO_PATH2CLEO=${path2CLEO}

if [ "${path2CLEO}" == "" ]; then
  echo "Please provide path to CLEO source directory"
  exit 1
fi

experiments_script="${path2CLEO}/scripts_2/common/experiments.sh"
if [ ! -f "${experiments_script}" ]; then
  echo "Error: experiments script not found at ${experiments_script}"
  exit 1
fi
source "${experiments_script}"
load_experiment_config "$5" "$6" "${experiment}"

yacyaxtroot=${7:-${HOME}/yacyaxt/${compilername}}
enabledebug=${8:-false}
make_clean=${9:-false}
stacksize_limit=${10:-204800}

if [[ "${buildtype}" == "cuda" && "${compilername}" != "gcc" ]]; then
  echo "Error: CUDA build on Levante requires compilername='gcc'."
  exit 1
fi
### ---------------------------------------------------- ###

### ----------------- export inputs -------------------- ###
export CLEO_YACYAXTROOT=${yacyaxtroot}
export CLEO_ENABLEDEBUG=${enabledebug}
export CLEO_MAKE_JOBS=32
export CLEO_PYTHON=${CLEO_PYTHON:-${path2CLEO}/.venv/bin/python3}

### ---------------------------------------------------- ###

### -------------------- check inputs ------------------ ###
source "${path2CLEO}/scripts_2/common/check_inputs.sh"

check_args_not_empty "${CLEO_BUILDTYPE}" "${CLEO_COMPILERNAME}" "${CLEO_PATH2CLEO}" \
                     "${CLEO_PATH2BUILD}" "${CLEO_BUILD_FLAGS}" "${CLEO_YACYAXTROOT}" \
                     "${CLEO_ENABLEDEBUG}" "${stacksize_limit}"
### ---------------------------------------------------- ###

print_config_script="${path2CLEO}/scripts_2/common/print_configuration.sh"
if [ -f "${print_config_script}" ]; then
  source "${print_config_script}"
  print_configuration "${experiment}"
fi

### --------------------- build CLEO ------------------- ###
buildcmd="${CLEO_PATH2CLEO}/scripts_2/common/build_cleo.sh"
if [ ! -f "${buildcmd}" ]; then
  echo "Error: build script not found at ${buildcmd}"
  exit 1
fi
echo "${buildcmd}"
source "${buildcmd}"
build_cleo
### ---------------------------------------------------- ###

### ---------------- compile experiment ---------------- ###
compilecmd="${CLEO_PATH2CLEO}/scripts_2/common/compile_cleo.sh"
echo "${compilecmd} \"${executables}\" ${make_clean}"
source "${compilecmd}"
compile_cleo "${executables}" "${make_clean}"
### ---------------------------------------------------- ###

### ------------- load runtime environment ------------- ###
runtime_settings_script="${path2CLEO}/scripts_2/levante/runtime_settings.sh"
if [ ! -f "${runtime_settings_script}" ]; then
  echo "Error: runtime settings script not found at ${runtime_settings_script}"
  exit 1
fi
source "${runtime_settings_script}"
configure_machine_runtime_settings "${stacksize_limit}"
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
