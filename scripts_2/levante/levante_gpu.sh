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

${SLURM_SUBMIT_DIR}/scripts_2/levante/build_compile_run_plot_cleo.sh constthermo2d cuda
