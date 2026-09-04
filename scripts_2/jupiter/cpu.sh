#!/bin/bash

#SBATCH --job-name=cleo_cpu
#SBATCH --partition=booster
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=288
#SBATCH --time=00:30:00
#SBATCH --account=xspies
#SBATCH --output=./cleo_cpu.%j.out
#SBATCH --error=./cleo_cpu.%j.out

set -e

# ============================================================
# Environment
# ============================================================

source /etc/profile

# Required if the .venv does not use the default environment.
module load Stages/2026 Python/3.13.5

export CLEO_PATH2CLEO="${SLURM_SUBMIT_DIR:-$(pwd)}"

# Use the CLEO uv environment.
export CLEO_PYTHON="${CLEO_PYTHON:-${CLEO_PATH2CLEO}/.venv/bin/python3}"

# YAC/YAXT installation.
export CLEO_YACYAXTROOT="${CLEO_YACYAXTROOT:-${HOME}/yacyaxt/gcc}"

# ============================================================
# Configuration
# ============================================================

mode="${1:-run}"

experiment="constthermo2d"
buildtype="openmp"
compilername="gcc"

# ============================================================
# Validation
# ============================================================

if [[ ! -x "${CLEO_PYTHON}" ]]; then
    echo "Error: CLEO Python executable not found:"
    echo "  ${CLEO_PYTHON}"
    exit 1
fi

source "${CLEO_PATH2CLEO}/scripts_2/common/check_inputs.sh"

check_args_not_empty \
    "${CLEO_PATH2CLEO}" \
    "${CLEO_PYTHON}" \
    "${CLEO_YACYAXTROOT}"

# ============================================================
# Helper functions
# ============================================================

run_cleo() {
    echo
    echo "=== Running ${experiment} (${buildtype}, ${compilername}) ==="
    echo

    srun --exclusive \
        --ntasks=1 \
        --cpus-per-task="${SLURM_CPUS_PER_TASK}" \
        "${CLEO_PATH2CLEO}/scripts_2/jupiter/build_compile_run_plot_cleo.sh" \
        "${experiment}" \
        "${buildtype}" \
        "${compilername}" \
        "${CLEO_PATH2CLEO}"
}

build_cleo() {
    echo
    echo "=== Building ${experiment} ==="
    echo

    "${CLEO_PATH2CLEO}/scripts_2/jupiter/build_compile_run_plot_cleo.sh" \
        "${experiment}" \
        "${buildtype}" \
        "${compilername}" \
        "${CLEO_PATH2CLEO}" \
        "" \
        "" \
        "${CLEO_YACYAXTROOT}" \
        false \
        false \
        204800 \
        build,compile
}

# ============================================================
# Main
# ============================================================

case "${mode}" in

    build)
        build_cleo
        ;;

    run)
        run_cleo
        ;;

    *)
        echo "Usage: $0 [build|run]"
        exit 1
        ;;

esac
