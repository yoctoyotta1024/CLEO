#!/bin/bash

# Slurm resources for the single-GPU CLEO job.
#SBATCH --job-name=cleo_gpu
#SBATCH --partition=booster
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=72
#SBATCH --gpus-per-node=1
#SBATCH --time=00:30:00
#SBATCH --account=xspies
#SBATCH --output=./cleo_gpu.%j.out
#SBATCH --error=./cleo_gpu.%j.out

set -e

# Load the site environment and the Python runtime used by the helper script.
source /etc/profile
module load Stages/2026 Python/3.13.5

# Resolve paths from the submission directory while allowing site-specific overrides.
export CLEO_PATH2CLEO="${SLURM_SUBMIT_DIR:-$(pwd)}"
export CLEO_PYTHON="${CLEO_PYTHON:-${CLEO_PATH2CLEO}/.venv/bin/python3}"
export CLEO_YACYAXTROOT="${CLEO_YACYAXTROOT:-${HOME}/yacyaxt/gcc}"

# Use run mode unless the caller explicitly requests a build.
mode="${1:-run}"
experiment="constthermo2d"
buildtype="cuda"
compilername="gcc"

source "${CLEO_PATH2CLEO}/scripts_2/common/check_inputs.sh"
check_args_not_empty "${CLEO_PATH2CLEO}" "${CLEO_PYTHON}" "${CLEO_YACYAXTROOT}"

# The build mode performs only the build and compile stages.
# The run mode omits the stages argument, so the helper defaults to all:
# build, compile, run, and plot.
case "${mode}" in
    build)
        "${CLEO_PATH2CLEO}/scripts_2/jupiter/build_compile_run_plot_cleo.sh" \
            "${experiment}" "${buildtype}" "${compilername}" "${CLEO_PATH2CLEO}" \
            "" "" "${CLEO_YACYAXTROOT}" false false 204800 build,compile
        ;;
    run)
        echo "=== Running ${experiment} (${buildtype}, ${compilername}) ==="
        "${CLEO_PATH2CLEO}/scripts_2/jupiter/build_compile_run_plot_cleo.sh" \
            "${experiment}" "${buildtype}" "${compilername}" "${CLEO_PATH2CLEO}"
        ;;
    *)
        echo "Usage: $0 [build|run]"
        exit 1
        ;;
esac
