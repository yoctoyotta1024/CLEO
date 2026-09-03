#!/bin/bash

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

source /etc/profile
module load Stages/2026 Python/3.13.5

export CLEO_PATH2CLEO="${SLURM_SUBMIT_DIR:-$(pwd)}"
export CLEO_PYTHON="${CLEO_PYTHON:-${CLEO_PATH2CLEO}/.venv/bin/python3}"
export CLEO_YACYAXTROOT="${CLEO_YACYAXTROOT:-${HOME}/yacyaxt/gcc}"

mode="${1:-run}"
experiment="constthermo2d"
buildtype="cuda"
compilername="gcc"

if [[ ! -x "${CLEO_PYTHON}" ]]; then
    echo "Error: CLEO Python executable not found: ${CLEO_PYTHON}"
    exit 1
fi

source "${CLEO_PATH2CLEO}/scripts_2/common/check_inputs.sh"
check_args_not_empty "${CLEO_PATH2CLEO}" "${CLEO_PYTHON}" "${CLEO_YACYAXTROOT}"

run_cleo() {
    echo "=== Running ${experiment} (${buildtype}, ${compilername}) ==="
    srun --exclusive \
        --ntasks=1 \
        --cpus-per-task="${SLURM_CPUS_PER_TASK:-72}" \
        --gpus-per-task=1 \
        "${CLEO_PATH2CLEO}/scripts_2/jupyter/build_compile_run_plot_cleo.sh" \
        "${experiment}" "${buildtype}" "${compilername}" "${CLEO_PATH2CLEO}"
}

case "${mode}" in
    build)
        "${CLEO_PATH2CLEO}/scripts_2/jupyter/build_compile_run_plot_cleo.sh" \
            "${experiment}" "${buildtype}" "${compilername}" "${CLEO_PATH2CLEO}" \
            "" "" "${CLEO_YACYAXTROOT}" false false 204800 build,compile
        ;;
    run)
        run_cleo
        ;;
    *)
        echo "Usage: $0 [build|run]"
        exit 1
        ;;
esac
