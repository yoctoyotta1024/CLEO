#!/bin/bash
#SBATCH --job-name=cleo_cpu
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=00:30:00
#SBATCH --mail-user=<YOUR_EMAIL>
#SBATCH --mail-type=FAIL
#SBATCH --account=<YOUR_ACCOUNT>
#SBATCH --output=./cleo_cpu.%j.out
#SBATCH --error=./cleo_cpu.%j.out

set -e
source /etc/profile
module purge
spack unload --all

CLEO_PATH2CLEO="${SLURM_SUBMIT_DIR:-$(pwd)}"

source "${CLEO_PATH2CLEO}/scripts_2/common/check_inputs.sh"
check_args_not_empty "${CLEO_PYTHON}" "${CLEO_YACYAXTROOT}"

run_experiment() {
  local experiment="$1"
  local buildtype="$2"
  local compilername="$3"

  echo "=== Running ${experiment} (${buildtype}, ${compilername}) ==="
  "${CLEO_PATH2CLEO}/scripts_2/levante/build_compile_run_plot_cleo.sh" \
    "${experiment}" "${buildtype}" "${compilername}" "${CLEO_PATH2CLEO}"
}

for entry in \
  "as2017 serial gcc" \
  "breakup serial gcc" \
  "constthermo2d openmp gcc" \
  "cuspbifurc threads gcc" \
  "divfree2d openmp gcc" \
  "eurec4a1d threads gcc" \
  "rainshaft1d threads gcc" \
  "shima2009 openmp gcc" \
  "python_bindings openmp gcc"; do
  set -- $entry
  run_experiment "$1" "$2" "$3"
done

# For a later MPI-heavy test, you can swap to:
# run_experiment fromfile openmp gcc
