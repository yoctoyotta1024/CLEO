#!/bin/bash
#SBATCH --job-name=cleo_gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --gpus=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=10G
#SBATCH --time=00:30:00
#SBATCH --mail-user=
#SBATCH --mail-type=FAIL
#SBATCH --account=
#SBATCH --output=./cleo_gpu.%j.out
#SBATCH --error=./cleo_gpu.%j.out

set -e
source /etc/profile
module purge
spack unload --all

repo_dir="${SLURM_SUBMIT_DIR:-$(pwd)}"

if [[ ! -f "${repo_dir}/scripts_2/levante/build_compile_run_plot_cleo.sh" ]]; then
  echo "Error: expected CLEO repo at '${repo_dir}' (missing scripts_2/levante/build_compile_run_plot_cleo.sh)."
  exit 1
fi

source "${repo_dir}/scripts_2/common/check_inputs.sh"
check_args_not_empty "${CLEO_PYTHON}" "${CLEO_YACYAXTROOT}"

"${repo_dir}/scripts_2/levante/build_compile_run_plot_cleo.sh" constthermo2d cuda gcc "${repo_dir}"
