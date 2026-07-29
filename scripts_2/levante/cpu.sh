#!/bin/bash
#SBATCH --job-name=cleo_cpu
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=32
#SBATCH --mem=10G
#SBATCH --time=00:30:00
#SBATCH --mail-user=
#SBATCH --mail-type=FAIL
#SBATCH --account=
#SBATCH --output=./cleo_cpu.%j.out
#SBATCH --error=./cleo_cpu.%j.out

set -e
source /etc/profile
module purge
spack unload --all

${SLURM_SUBMIT_DIR}/scripts_2/levante/build_compile_run_plot_cleo.sh fromfile openmp
