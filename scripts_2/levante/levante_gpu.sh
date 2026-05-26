#!/bin/bash
#SBATCH --job-name=build_compile_run_cleo
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
#SBATCH --output=./build_compile_run_cleo_out.%j.out
#SBATCH --error=./build_compile_run_cleo_err.%j.out

set -e
source /etc/profile
module purge
spack unload --all

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
${SCRIPT_DIR}/build_compile_run_plot_cleo.sh "$@"
