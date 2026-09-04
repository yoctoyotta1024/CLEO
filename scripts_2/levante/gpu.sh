#!/bin/bash
#SBATCH --job-name=cleo_gpu
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --gpus=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=00:30:00
#SBATCH --mail-user=
#SBATCH --mail-type=FAIL
#SBATCH --account=
#SBATCH --output=./cleo_gpu.%j.out
#SBATCH --error=./cleo_gpu.%j.out

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
  "as2017 cuda gcc" \
  "breakup cuda gcc" \
  "constthermo2d cuda gcc" \
  "cuspbifurc cuda gcc" \
  "divfree2d cuda gcc" \
  "eurec4a1d cuda gcc" \
  "rainshaft1d cuda gcc" \
  "shima2009 cuda gcc" \
  "python_bindings cuda gcc"; do
  set -- $entry
  run_experiment "$1" "$2" "$3"
done

# For a later MPI-heavy test, you can swap to:
# run_experiment fromfile openmp gcc

#"${CLEO_PATH2CLEO}/scripts_2/levante/build_compile_run_plot_cleo.sh" constthermo2d cuda gcc "${CLEO_PATH2CLEO}"
