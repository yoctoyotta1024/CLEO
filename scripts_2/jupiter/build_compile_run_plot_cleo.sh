#!/bin/bash

### ============================================================ ###
###                     jupiter script                          ###
### ============================================================ ###
###
### Usage:
###   ./build_compile_run_plot_cleo.sh [experiment] [buildtype] [compilername] \
###                                    [path2CLEO] [path2build] [build_flags] \
###                                    [yacyaxtroot] [enabledebug] [make_clean] \
###                                    [stacksize_limit] [steps]
###
### Arguments:
###   $1  experiment       Name of experiment                  (default: as2017)
###   $2  buildtype        serial | threads | openmp | cuda    (default: openmp)
###   $3  compilername     gcc                                (default: gcc)
###   $4  path2CLEO        Absolute path to CLEO source        (default: $HOME/CLEO)
###   $5  path2build       Absolute build path                 (default: experiment)
###   $6  build_flags      Extra CMake flags                  (default: experiment)
###   $7  yacyaxtroot      Path to YAC+YAXT installation       (default: $HOME/yacyaxt/gcc)
###   $8  enabledebug      true | false                       (default: false)
###   $9  make_clean        true | false                       (default: false)
###   $10 stacksize_limit  ulimit -s value (kB)               (default: 204800)
###   $11 steps             build,compile,run,plot,all         (default: all)
###
###   The run stage generates input files and runs the executable.
###   Set NSYS_PREFIX to profile only the executable run.
###
### Available examples (see common/experiments.sh for full details):
###   as2017  cuspbifurc  breakup  shima2009  constthermo2d  divfree2d
###   eurec4a1d  rainshaft1d  python_bindings
###   fromfile  fromfile_irreg  bubble3d
### ============================================================ ###

set -e

export CLEO_MACHINE=jupiter

experiment=${1:-as2017}
buildtype=${2:-openmp}
compilername=${3:-gcc}
path2CLEO=${4:-${HOME}/CLEO}
path2build=${5:-}
build_flags=${6:-}
yacyaxtroot=${7:-${HOME}/yacyaxt/gcc}
enabledebug=${8:-false}
make_clean=${9:-false}
stacksize_limit=${10:-204800}
steps=${11:-all}

if [[ ! -d "${path2CLEO}" ]]; then
    echo "Error: CLEO source directory not found: ${path2CLEO}"
    exit 1
fi
if [[ "${compilername}" != gcc ]]; then
    echo "Error: jupiter currently supports compilername='gcc'."
    exit 1
fi
if [[ "${steps}" != all ]]; then
    IFS=',' read -r -a requested_steps <<< "${steps}"
    for step in "${requested_steps[@]}"; do
        case "${step}" in
            build|compile|run|plot) ;;
            *)
                echo "Error: invalid step '${step}'. Use build, compile, run, plot, or all."
                exit 1
                ;;
        esac
    done
fi

export CLEO_BUILDTYPE=${buildtype}
export CLEO_COMPILERNAME=${compilername}
export CLEO_PATH2CLEO=${path2CLEO}
export CLEO_PATH2BUILD=${path2build}
export CLEO_BUILD_FLAGS=${build_flags}
export CLEO_YACYAXTROOT=${yacyaxtroot}
export CLEO_ENABLEDEBUG=${enabledebug}
export CLEO_MAKE_JOBS=${CLEO_MAKE_JOBS:-32}
export CLEO_PYTHON=${CLEO_PYTHON:-${path2CLEO}/.venv/bin/python3}

source "${path2CLEO}/scripts_2/common/experiments.sh"
load_experiment_config "${path2build}" "${build_flags}" "${experiment}"

step_enabled() {
    [[ "${steps}" == all || ",${steps}," == *",$1,"* ]]
}

source "${path2CLEO}/scripts_2/common/check_inputs.sh"
check_args_not_empty "${CLEO_BUILDTYPE}" "${CLEO_COMPILERNAME}" "${CLEO_PATH2CLEO}" \
                                         "${CLEO_PATH2BUILD}" "${CLEO_BUILD_FLAGS}" "${CLEO_YACYAXTROOT}" \
                                         "${CLEO_ENABLEDEBUG}" "${stacksize_limit}"

if [[ -f "${path2CLEO}/scripts_2/common/print_configuration.sh" ]]; then
    source "${path2CLEO}/scripts_2/common/print_configuration.sh"
    print_configuration "${experiment}"
fi

if step_enabled build; then
    source "${path2CLEO}/scripts_2/common/build_cleo.sh"
    build_cleo
fi

if step_enabled compile; then
    source "${path2CLEO}/scripts_2/common/compile_cleo.sh"
    compile_cleo "${executables}" "${make_clean}"
fi

if step_enabled run || step_enabled plot; then
    source "${path2CLEO}/scripts_2/jupiter/runtime_settings.sh"
    configure_machine_runtime_settings "${stacksize_limit}"

    if [[ ! -f "${pythonscript}" ]]; then
        echo "Error: Python script not found: ${pythonscript}"
        exit 1
    fi

    read -r -a experiment_args <<< "${script_args:-}"
    python_args=()
    for arg in "${experiment_args[@]}"; do
        case "${arg}" in
            --do_inputfiles|--do_run_executable|--do_plot_results) ;;
            *) python_args+=("${arg}") ;;
        esac
    done
fi

run_python_stage() {
    local flag="$1"
    local prefix=()

    if [[ "${flag}" == --do_run_executable && -n "${NSYS_PREFIX:-}" ]]; then
        read -r -a prefix <<< "${NSYS_PREFIX}"
        echo "Running (profiled): ${NSYS_PREFIX} ${CLEO_PYTHON} ${pythonscript}"
    else
        echo "Running: ${CLEO_PYTHON} ${pythonscript} ${flag}"
    fi

    "${prefix[@]}" "${CLEO_PYTHON}" "${pythonscript}" \
        "${path2CLEO}" "${CLEO_PATH2BUILD}" "${python_args[@]}" "${flag}"
}

if step_enabled run; then
    run_python_stage --do_inputfiles
    run_python_stage --do_run_executable
fi

if step_enabled plot; then
    run_python_stage --do_plot_results
fi
