#!/bin/bash

### ============================================================ ###
###        Shared build / compile / run / plot pipeline          ###
### ============================================================ ###
###
### Not run directly: each <machine>/build_compile_run_plot_cleo.sh sets
### the machine configuration below and then calls
###   build_compile_run_plot_cleo "$@"
###
### Machine configuration (set by the machine script):
###   CLEO_MACHINE                 machine name (folder in scripts/)
###   machine_default_buildtype    default for $2
###   machine_buildtypes=(...)     supported build types
###   machine_compilers=(...)      supported compiler names
###   machine_examples=(...)    supported examples
###   machine_default_stacksize    default for $10 ("" = leave unchanged)
###   machine_default_make_jobs    default for CLEO_MAKE_JOBS
###   machine_check_inputs()       optional extra validation
###
### Arguments (identical on every machine, all optional):
###   $1  example       Name of example                  (default: as2017)
###   $2  buildtype        see machine script                  (default: machine)
###   $3  compilername     see machine script                  (default: gcc)
###   $4  path2CLEO        Absolute path to CLEO source        (default: $CLEO_PATH2CLEO, else $HOME/CLEO)
###   $5  path2build       Exact build directory for the example
###                        (default: <path2CLEO>/build_xxx from examples.sh)
###   $6  build_flags      Extra CMake flags                   (default: example)
###   $7  yacyaxtroot      YAC+YAXT root, with yac and yaxt in <yacyaxtroot>/<compilername>
###                        (default: $CLEO_YACYAXTROOT, else $HOME/yacyaxt)
###   $8  enabledebug      true | false                        (default: $CLEO_ENABLEDEBUG, else false)
###   $9  stacksize_limit  ulimit -s value (kB)                (default: $CLEO_STACKSIZE_LIMIT, else machine)
###   $10 steps            configure,compile,inputs,run,plot,all (default: all)
###
### Environment:
###   CLEO_PYTHON     python to run the example scripts (default: <path2CLEO>/.venv/bin/python3)
###   CLEO_MAKE_JOBS  parallel make jobs                   (default: machine)
### ============================================================ ###

set -e

step_enabled() {
  [[ "${steps}" == all || ",${steps}," == *",$1,"* ]]
}

run_python_stage() {
  local flag="$1"

  echo "Running: ${CLEO_PYTHON} ${pythonscript} ${flag}"
  "${CLEO_PYTHON}" "${pythonscript}" \
    "${path2CLEO}" "${CLEO_PATH2BUILD}" "${python_args[@]}" "${flag}"
}

build_compile_run_plot_cleo() {
  local common_dir
  common_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
  local machine_dir="${common_dir}/../${CLEO_MACHINE}"

  source "${common_dir}/check_inputs.sh"
  check_machine

  ### ---------------- read arguments ------------------ ###
  example=${1:-as2017}
  buildtype=${2:-${machine_default_buildtype}}
  compilername=${3:-gcc}
  path2CLEO=${4:-${CLEO_PATH2CLEO:-${HOME}/CLEO}}
  local path2build_override=${5:-}
  local build_flags_override=${6:-}
  # yac and yaxt for each compiler are in <yacyaxtroot>/<compilername>
  yacyaxtroot=${7:-${CLEO_YACYAXTROOT:-${HOME}/yacyaxt}}
  yacyaxtroot="${yacyaxtroot%/}/${compilername}"
  enabledebug=${8:-${CLEO_ENABLEDEBUG:-false}}
  stacksize_limit=${9:-${CLEO_STACKSIZE_LIMIT:-${machine_default_stacksize}}}
  steps=${10:-all}
  ### ---------------------------------------------------- ###

  ### ------------------ check arguments --------------- ###
  if [[ ! -d "${path2CLEO}" ]]; then
    echo "Error: CLEO source directory not found: ${path2CLEO}"
    exit 1
  fi
  check_value_in_list example "${example}" "${machine_examples[@]}"
  check_value_in_list buildtype "${buildtype}" "${machine_buildtypes[@]}"
  check_value_in_list compilername "${compilername}" "${machine_compilers[@]}"
  check_value_in_list enabledebug "${enabledebug}" true false
  check_steps "${steps}"
  if declare -F machine_check_inputs >/dev/null; then
    machine_check_inputs
  fi
  ### ---------------------------------------------------- ###

  ### ----------------- export inputs ------------------- ###
  export CLEO_BUILDTYPE=${buildtype}
  export CLEO_COMPILERNAME=${compilername}
  export CLEO_PATH2CLEO=${path2CLEO}
  export CLEO_YACYAXTROOT=${yacyaxtroot}
  export CLEO_ENABLEDEBUG=${enabledebug}
  export CLEO_MAKE_JOBS=${CLEO_MAKE_JOBS:-${machine_default_make_jobs}}
  export CLEO_PYTHON=${CLEO_PYTHON:-${path2CLEO}/.venv/bin/python3}

  source "${common_dir}/examples.sh"
  load_example_config "${path2build_override}" "${build_flags_override}" "${example}"
  check_source_and_build_paths

  check_args_not_empty "${CLEO_BUILDTYPE}" "${CLEO_COMPILERNAME}" "${CLEO_PATH2CLEO}" \
                       "${CLEO_PATH2BUILD}" "${CLEO_BUILD_FLAGS}" "${CLEO_YACYAXTROOT}" \
                       "${CLEO_ENABLEDEBUG}"

  source "${common_dir}/print_configuration.sh"
  print_configuration "${example}"
  ### ---------------------------------------------------- ###

  ### --------------------- build CLEO ------------------ ###
  if step_enabled configure; then
    source "${common_dir}/build_cleo.sh"
    build_cleo
  fi
  ### ---------------------------------------------------- ###

  ### ---------------- compile example -------------- ###
  if step_enabled compile; then
    if ! step_enabled configure; then
      # the configure step loads the machine's toolchain and YAC's dependencies
      # (modules), so without it load them here before compiling; the linker
      # needs YAC's dependencies (e.g. NetCDF -> HDF5, PnetCDF) too
      source "${machine_dir}/build_flags.sh"
      configure_machine_build_flags
      source "${machine_dir}/helpers/build_yac.sh"
      configure_machine_yac_flags
    fi
    source "${common_dir}/compile_cleo.sh"
    compile_cleo "${executables}"
  fi
  ### ---------------------------------------------------- ###

  ### ------ load runtime environment + python args ----- ###
  if step_enabled inputs || step_enabled run || step_enabled plot; then
    source "${machine_dir}/runtime_settings.sh"
    configure_machine_runtime_settings "${stacksize_limit}"

    if [[ ! -f "${pythonscript}" ]]; then
      echo "Error: Python script not found: ${pythonscript}"
      exit 1
    fi
    if [[ ! -x "${CLEO_PYTHON}" ]] && ! command -v "${CLEO_PYTHON}" &>/dev/null; then
      echo "Error: CLEO_PYTHON not found: ${CLEO_PYTHON}"
      exit 1
    fi

    # example's python args (from examples.sh); each stage adds its own --do_* flag
    python_args=()
    read -r -a python_args <<< "${script_args:-}"
  fi
  ### ---------------------------------------------------- ###

  ### --------- inputs / run / plot example --------- ###
  if step_enabled inputs; then
    run_python_stage --do_inputfiles
  fi

  if step_enabled run; then
    run_python_stage --do_run_executable
  fi

  if step_enabled plot; then
    run_python_stage --do_plot_results
  fi
  ### ---------------------------------------------------- ###
}
