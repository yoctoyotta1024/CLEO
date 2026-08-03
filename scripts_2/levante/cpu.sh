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

repo_dir="${SLURM_SUBMIT_DIR:-$(pwd)}"

if [[ ! -f "${repo_dir}/scripts_2/levante/build_compile_run_plot_cleo.sh" ]]; then
  echo "Error: expected CLEO repo at '${repo_dir}' (missing scripts_2/levante/build_compile_run_plot_cleo.sh)."
  exit 1
fi

source "${repo_dir}/scripts_2/common/check_inputs.sh"
check_args_not_empty "${CLEO_PYTHON}" "${CLEO_YACYAXTROOT}"

# kokkostools is part of this matrix; require explicit path or opt-in auto-install inputs.
if [[ -z "${CLEO_KOKKOSTOOLS}" && "${CLEO_AUTO_INSTALL_KOKKOSTOOLS}" != "true" ]]; then
  echo "Error: set CLEO_KOKKOSTOOLS, or set CLEO_AUTO_INSTALL_KOKKOSTOOLS=true."
  exit 1
fi
if [[ -z "${CLEO_KOKKOSTOOLS}" && "${CLEO_AUTO_INSTALL_KOKKOSTOOLS}" == "true" && -z "${CLEO_KOKKOSTOOLS_REPO_PARENT}" ]]; then
  echo "Error: CLEO_AUTO_INSTALL_KOKKOSTOOLS=true requires CLEO_KOKKOSTOOLS_REPO_PARENT."
  exit 1
fi

run_experiment() {
  local experiment="$1"
  local buildtype="$2"
  local compilername="$3"

  echo "=== Running ${experiment} (${buildtype}, ${compilername}) ==="
  "${repo_dir}/scripts_2/levante/build_compile_run_plot_cleo.sh" \
    "${experiment}" "${buildtype}" "${compilername}" "${repo_dir}"
}

for entry in \
  "as2017 serial gcc" \
  "breakup serial gcc" \
  "constthermo2d serial gcc" \
  "cuspbifurc serial gcc" \
  "divfree2d serial gcc" \
  "eurec4a1d serial gcc" \
  "kokkostools serial gcc" \
  "python_bindings serial gcc" \
  "rainshaft1d serial gcc" \
  "shima2009 serial gcc" \
  "as2017 openmp gcc" \
  "breakup openmp gcc" \
  "divfree2d openmp gcc" \
  "rainshaft1d openmp gcc"; do
  set -- $entry
  run_experiment "$1" "$2" "$3"
done

# For a later MPI-heavy test, you can swap to:
# run_experiment fromfile openmp gcc
