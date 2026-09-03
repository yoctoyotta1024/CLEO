```bash
#!/bin/bash

### ============================================================ ###
###                     Jupyter script                          ###
### ============================================================ ###
#
# Usage:
#   ./build_compile_run_plot_cleo.sh [experiment] [buildtype] [compilername] \
#                                    [path2CLEO] [path2build] [build_flags] \
#                                    [yacyaxtroot] [enabledebug] [make_clean] \
#                                    [stacksize_limit] [steps]
#
# All arguments are optional. Defaults come from the environment
# and experiments.sh.
#
# steps is a comma-separated list containing:
#   build, compile, run, plot, all
#
# Default:
#   all
#
# Examples:
#   as2017  cuspbifurc  breakup  shima2009  constthermo2d  divfree2d
#   eurec4a1d  rainshaft1d  python_bindings
#   fromfile  fromfile_irreg  bubble3d
#
### ============================================================ ###

set -e

# ------------------------------------------------------------------
# Machine configuration
# ------------------------------------------------------------------

export CLEO_MACHINE="jupyter"

# ------------------------------------------------------------------
# Arguments
# ------------------------------------------------------------------

experiment=${1:-as2017}
buildtype=${2:-openmp}
compilername=${3:-gcc}
path2CLEO=${4:-${HOME}/CLEO}
path2build=${5:-}
build_flags=${6:-}
yacyaxtroot=${7:-${HOME}/yacyaxt/${compilername}}
enabledebug=${8:-false}
make_clean=${9:-false}
stacksize_limit=${10:-204800}
steps=${11:-all}

# ------------------------------------------------------------------
# Basic validation
# ------------------------------------------------------------------

if [[ ! -d "${path2CLEO}" ]]; then
    echo "Error: CLEO source directory not found:"
    echo "  ${path2CLEO}"
    exit 1
fi

if [[ "${compilername}" != "gcc" ]]; then
    echo "Error: Jupyter currently supports compilername='gcc'."
    exit 1
fi

case ",${steps}," in
    *,build,*|*,compile,*|*,run,*|*,plot,*|,all,)
        ;;
    *)
        echo "Error: invalid steps '${steps}'."
        echo "       Use: build, compile, run, plot, or all."
        exit 1
        ;;
esac

# ------------------------------------------------------------------
# CLEO environment
# ------------------------------------------------------------------

export CLEO_BUILDTYPE="${buildtype}"
export CLEO_COMPILERNAME="${compilername}"
export CLEO_PATH2CLEO="${path2CLEO}"
export CLEO_PATH2BUILD="${path2build}"
export CLEO_BUILD_FLAGS="${build_flags}"
export CLEO_YACYAXTROOT="${yacyaxtroot}"
export CLEO_ENABLEDEBUG="${enabledebug}"
export CLEO_MAKE_JOBS="${CLEO_MAKE_JOBS:-32}"

# Prefer the CLEO uv environment.
if [[ -z "${CLEO_PYTHON:-}" ]]; then
    CLEO_PYTHON="${path2CLEO}/.venv/bin/python3"
fi

export CLEO_PYTHON

if [[ ! -x "${CLEO_PYTHON}" ]]; then
    echo "Error: Python executable not found:"
    echo "  ${CLEO_PYTHON}"
    echo
    echo "Create the CLEO environment first, e.g.:"
    echo "  uv sync"
    exit 1
fi

# ------------------------------------------------------------------
# Load experiment configuration
# ------------------------------------------------------------------

experiments_script="${path2CLEO}/scripts_2/common/experiments.sh"

if [[ ! -f "${experiments_script}" ]]; then
    echo "Error: experiments script not found:"
    echo "  ${experiments_script}"
    exit 1
fi

source "${experiments_script}"

load_experiment_config "${path2build}" "${build_flags}" "${experiment}"

# ------------------------------------------------------------------
# Check configuration
# ------------------------------------------------------------------

source "${path2CLEO}/scripts_2/common/check_inputs.sh"

check_args_not_empty \
    "${CLEO_BUILDTYPE}" \
    "${CLEO_COMPILERNAME}" \
    "${CLEO_PATH2CLEO}" \
    "${CLEO_PATH2BUILD}" \
    "${CLEO_BUILD_FLAGS}" \
    "${CLEO_YACYAXTROOT}" \
    "${CLEO_ENABLEDEBUG}" \
    "${stacksize_limit}"

# ------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------

step_enabled() {
    [[ "${steps}" == "all" || ",${steps}," == *",$1,"* ]]
}

run_python_script() {
    local flag="$1"
    shift

    echo
    echo "Running:"
    echo "  ${CLEO_PYTHON} ${pythonscript} ${path2CLEO} ${CLEO_PATH2BUILD} $* ${flag}"
    echo

    "${CLEO_PYTHON}" \
        "${pythonscript}" \
        "${path2CLEO}" \
        "${CLEO_PATH2BUILD}" \
        "${filtered_args[@]}" \
        "${flag}"
}

# ------------------------------------------------------------------
# Print configuration
# ------------------------------------------------------------------

print_config_script="${path2CLEO}/scripts_2/common/print_configuration.sh"

if [[ -f "${print_config_script}" ]]; then
    source "${print_config_script}"
    print_configuration "${experiment}"
fi

# ------------------------------------------------------------------
# Build
# ------------------------------------------------------------------

if step_enabled build; then

    build_script="${path2CLEO}/scripts_2/common/build_cleo.sh"

    if [[ ! -f "${build_script}" ]]; then
        echo "Error: build script not found:"
        echo "  ${build_script}"
        exit 1
    fi

    source "${build_script}"
    build_cleo
fi

# ------------------------------------------------------------------
# Compile
# ------------------------------------------------------------------

if step_enabled compile; then

    compile_script="${path2CLEO}/scripts_2/common/compile_cleo.sh"

    if [[ ! -f "${compile_script}" ]]; then
        echo "Error: compile script not found:"
        echo "  ${compile_script}"
        exit 1
    fi

    source "${compile_script}"
    compile_cleo "${executables}" "${make_clean}"
fi

# ------------------------------------------------------------------
# Runtime configuration
# ------------------------------------------------------------------

if step_enabled run || step_enabled plot; then

    runtime_settings_script="${path2CLEO}/scripts_2/jupyter/runtime_settings.sh"

    if [[ ! -f "${runtime_settings_script}" ]]; then
        echo "Error: runtime settings script not found:"
        echo "  ${runtime_settings_script}"
        exit 1
    fi

    source "${runtime_settings_script}"
    configure_machine_runtime_settings "${stacksize_limit}"
fi

# ------------------------------------------------------------------
# Python script validation
# ------------------------------------------------------------------

if step_enabled run || step_enabled plot; then

    if [[ -z "${pythonscript:-}" ]]; then
        echo "Error: no Python script provided for experiment '${experiment}'."
        exit 1
    fi

    if [[ ! -f "${pythonscript}" ]]; then
        echo "Error: Python script not found:"
        echo "  ${pythonscript}"
        exit 1
    fi
fi

# ------------------------------------------------------------------
# Clean experiment arguments
#
# The Python driver may contain these flags already. Remove them here
# so that this script controls which stage is executed.
# ------------------------------------------------------------------

read -r -a base_args <<< "${script_args:-}"

filtered_args=()

for arg in "${base_args[@]}"; do
    case "${arg}" in
        --do_inputfiles|--do_run_executable|--do_plot_results)
            ;;
        *)
            filtered_args+=("${arg}")
            ;;
    esac
done

# ------------------------------------------------------------------
# Run
# ------------------------------------------------------------------

if step_enabled run; then

    # Generate input files.
    run_python_script "--do_inputfiles"

    # Run executable, optionally through Nsight Systems.
    echo
    echo "Running executable..."

    if [[ -n "${NSYS_PREFIX:-}" ]]; then

        read -r -a nsys_cmd <<< "${NSYS_PREFIX}"

        echo "Profiling with:"
        echo "  ${NSYS_PREFIX}"

        "${nsys_cmd[@]}" \
            "${CLEO_PYTHON}" \
            "${pythonscript}" \
            "${path2CLEO}" \
            "${CLEO_PATH2BUILD}" \
            "${filtered_args[@]}" \
            --do_run_executable

    else

        "${CLEO_PYTHON}" \
            "${pythonscript}" \
            "${path2CLEO}" \
            "${CLEO_PATH2BUILD}" \
            "${filtered_args[@]}" \
            --do_run_executable

    fi
fi

# ------------------------------------------------------------------
# Plot
# ------------------------------------------------------------------

if step_enabled plot; then
    run_python_script "--do_plot_results"
fi
```
