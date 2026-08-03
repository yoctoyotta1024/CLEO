#!/bin/bash

set -e

### ============================================================ ###
###                      Vanilla script                          ###
### ============================================================ ###
###
### Usage:
###   ./build_compile_run_plot_cleo.sh [experiment] [buildtype] [compilername] \
###                                    [path2CLEO] [path2build] [build_flags] \
###                                    [yacyaxtroot] [enabledebug] [make_clean]
###
###   All arguments are optional — defaults come from environment variables
###   and experiments.sh.
###
### Arguments:
###   $1  experiment     Name of experiment to run       (default: set below)
###   $2  buildtype      serial | threads | openmp       (default: from experiment)
###   $3  compilername   gcc                             (default: from experiment)
###   $4  path2CLEO      Absolute path to CLEO source    (default: $CLEO_PATH2CLEO)
###   $5  path2build     Absolute path for build dir     (default: build_<experiment>)
###   $6  build_flags    Extra CMake flags               (default: from experiment)
###   $7  yacyaxtroot    Path to YAC+YAXT install        (default: $CLEO_YACYAXTROOT)
###   $8  enabledebug    true | false                    (default: false)
###   $9  make_clean     true | false                    (default: false)
###
### Required env vars (unless passed as args):
###   CLEO_PATH2CLEO
###   CLEO_PYTHON
###   CLEO_YACYAXTROOT
###
### Available examples (see common/experiments.sh for full details):
###   as2017  breakup  constthermo2d  cuspbifurc  divfree2d
###   eurec4a1d  kokkostools  python_bindings  rainshaft1d  shima2009
###
###   Note: fromfile, fromfile_irreg, and bubble3d require Levante (SLURM/YAC)
### ============================================================ ###

#step 1: configure
export CLEO_MACHINE="vanilla"

experiment=${1:-"as2017"}
buildtype=${2:-"serial"}
compilername=${3:-gcc}
path2CLEO=${4:-${CLEO_PATH2CLEO}}

# Export build context early so experiment resolution can use it (e.g. kokkostools paths).
export CLEO_BUILDTYPE=${buildtype}
export CLEO_COMPILERNAME=${compilername}
export CLEO_PATH2CLEO=${path2CLEO}

### --- validate experiment is supported on vanilla --- ###
vanilla_experiments=(
  as2017 breakup constthermo2d cuspbifurc divfree2d
  eurec4a1d kokkostools python_bindings rainshaft1d shima2009
)
valid=false
for e in "${vanilla_experiments[@]}"; do
  [[ "${experiment}" == "${e}" ]] && valid=true && break
done
if [[ "${valid}" == false ]]; then
  echo "Error: experiment '${experiment}' is not supported on vanilla."
  echo "Supported experiments: ${vanilla_experiments[*]}"
  exit 1
fi
### --------------------------------------------------- ###

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

yacyaxtroot=${7:-${CLEO_YACYAXTROOT}}
enabledebug=${8:-false}
make_clean=${9:-false}

### ----------------- export inputs -------------------- ###
export CLEO_YACYAXTROOT=${yacyaxtroot}
export CLEO_ENABLEDEBUG=${enabledebug}
### ---------------------------------------------------- ###

### -------------------- check inputs ------------------ ###
source "${path2CLEO}/scripts_2/common/check_inputs.sh"

check_args_not_empty "${CLEO_BUILDTYPE}" "${CLEO_COMPILERNAME}" "${CLEO_PATH2CLEO}" \
                     "${CLEO_PATH2BUILD}" "${CLEO_BUILD_FLAGS}" "${CLEO_YACYAXTROOT}" \
                     "${CLEO_ENABLEDEBUG}"

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

### ---------------- compile experiment --------------- ###
compilecmd="${CLEO_PATH2CLEO}/scripts_2/common/compile_cleo.sh"
echo "${compilecmd} \"${executables}\" ${make_clean}"
source "${compilecmd}"
compile_cleo "${executables}" "${make_clean}"
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
