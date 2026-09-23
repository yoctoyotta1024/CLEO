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
###   CLEO_MACHINE                 machine name (folder in scripts_2/)
###   machine_default_buildtype    default for $2
###   machine_buildtypes=(...)     supported build types
###   machine_compilers=(...)      supported compiler names
###   machine_experiments=(...)    supported experiments
###   machine_default_stacksize    default for $10 ("" = leave unchanged)
###   machine_default_make_jobs    default for CLEO_MAKE_JOBS
###   machine_check_inputs()       optional extra validation
###
### Arguments (identical on every machine, all optional):
###   $1  experiment       Name of experiment                  (default: as2017)
###   $2  buildtype        see machine script                  (default: machine)
###   $3  compilername     see machine script                  (default: gcc)
###   $4  path2CLEO        Absolute path to CLEO source        (default: $CLEO_PATH2CLEO, else $HOME/CLEO)
###   $5  path2build       Build root folder; the experiment builds in
###                        <path2build>/build_xxx  (default: path2CLEO)
###   $6  build_flags      Extra CMake flags                   (default: experiment)
###   $7  yacyaxtroot      Path to YAC+YAXT installation       (default: $CLEO_YACYAXTROOT, else $HOME/yacyaxt/<compilername>)
###   $8  enabledebug      true | false                        (default: false)
###   $9  make_clean       true | false                        (default: false)
###   $10 stacksize_limit  ulimit -s value (kB)                (default: machine)
###   $11 steps            build,compile,run,plot,all          (default: all)
###
### Environment:
###   CLEO_PYTHON     python to run the experiment scripts (default: <path2CLEO>/.venv/bin/python3)
###   CLEO_MAKE_JOBS  parallel make jobs                   (default: machine)
###   NSYS_PREFIX     command prefix (e.g. a profiler) for the executable run only
### ============================================================ ###

set -e

step_enabled() {
  [[ "${steps}" == all || ",${steps}," == *",$1,"* ]]
}

run_python_stage() {
  local flag="$1"
  local prefix=()

  if [[ "${flag}" == --do_run_executable && -n "${NSYS_PREFIX:-}" ]]; then
    read -r -a prefix <<< "${NSYS_PREFIX}"
    echo "Running (profiled): ${NSYS_PREFIX} ${CLEO_PYTHON} ${pythonscript} ${flag}"
  else
    echo "Running: ${CLEO_PYTHON} ${pythonscript} ${flag}"
  fi

  "${prefix[@]}" "${CLEO_PYTHON}" "${pythonscript}" \
    "${path2CLEO}" "${CLEO_PATH2BUILD}" "${python_args[@]}" "${flag}"
}

build_compile_run_plot_cleo() {
  local common_dir
  common_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
  local machine_dir="${common_dir}/../${CLEO_MACHINE}"

  source "${common_dir}/check_inputs.sh"
  check_machine

  ### ---------------- read arguments ------------------ ###
  experiment=${1:-as2017}
  buildtype=${2:-${machine_default_buildtype}}
  compilername=${3:-gcc}
  path2CLEO=${4:-${CLEO_PATH2CLEO:-${HOME}/CLEO}}
  local path2build_override=${5:-}
  local build_flags_override=${6:-}
  yacyaxtroot=${7:-${CLEO_YACYAXTROOT:-${HOME}/yacyaxt/${compilername}}}
  enabledebug=${8:-false}
  make_clean=${9:-false}
  stacksize_limit=${10:-${machine_default_stacksize}}
  steps=${11:-all}
  ### ---------------------------------------------------- ###

  ### ------------------ check arguments --------------- ###
  if [[ ! -d "${path2CLEO}" ]]; then
    echo "Error: CLEO source directory not found: ${path2CLEO}"
    exit 1
  fi
  check_value_in_list experiment "${experiment}" "${machine_experiments[@]}"
  check_value_in_list buildtype "${buildtype}" "${machine_buildtypes[@]}"
  check_value_in_list compilername "${compilername}" "${machine_compilers[@]}"
  check_value_in_list enabledebug "${enabledebug}" true false
  check_value_in_list make_clean "${make_clean}" true false
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

  source "${common_dir}/experiments.sh"
  load_experiment_config "${path2build_override}" "${build_flags_override}" "${experiment}"

  check_args_not_empty "${CLEO_BUILDTYPE}" "${CLEO_COMPILERNAME}" "${CLEO_PATH2CLEO}" \
                       "${CLEO_PATH2BUILD}" "${CLEO_BUILD_FLAGS}" "${CLEO_YACYAXTROOT}" \
                       "${CLEO_ENABLEDEBUG}"

  source "${common_dir}/print_configuration.sh"
  print_configuration "${experiment}"
  ### ---------------------------------------------------- ###

  ### --------------------- build CLEO ------------------ ###
  if step_enabled build; then
    source "${common_dir}/build_cleo.sh"
    build_cleo
  fi
  ### ---------------------------------------------------- ###

  ### ---------------- compile experiment -------------- ###
  if step_enabled compile; then
    source "${common_dir}/compile_cleo.sh"
    compile_cleo "${executables}" "${make_clean}"
  fi
  ### ---------------------------------------------------- ###

  ### ------ load runtime environment + python args ----- ###
  if step_enabled run || step_enabled plot; then
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

    # strip stage flags from the experiment's args; each stage adds its own
    local experiment_args
    read -r -a experiment_args <<< "${script_args:-}"
    python_args=()
    local arg
    for arg in "${experiment_args[@]}"; do
      case "${arg}" in
        --do_inputfiles|--do_run_executable|--do_plot_results) ;;
        *) python_args+=("${arg}") ;;
      esac
    done
  fi
  ### ---------------------------------------------------- ###

  ### ------------- run / plot experiment -------------- ###
  if step_enabled run; then
    run_python_stage --do_inputfiles
    run_python_stage --do_run_executable
  fi

  if step_enabled plot; then
    run_python_stage --do_plot_results
  fi
  ### ---------------------------------------------------- ###
}
