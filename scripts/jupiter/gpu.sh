#!/bin/bash
#SBATCH --job-name=cleo_gpu
#SBATCH --partition=booster
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=72
#SBATCH --gpus-per-node=1
#SBATCH --time=00:30:00
#SBATCH --account=<YOUR_ACCOUNT>
#SBATCH --output=./cleo_gpu.%j.out
#SBATCH --error=./cleo_gpu.%j.out

### ============================================================ ###
###                    Jupiter GPU job script                    ###
### ============================================================ ###
###
### Usage:
###   ./scripts/jupiter/gpu.sh configure,compile [example] [buildtype] [compilername]
###   sbatch scripts/jupiter/gpu.sh inputs,run,plot [example] [buildtype] [compilername]
###   (buildtype is always cuda; any other value is ignored)
###
###   Run from (or submit from) the CLEO root directory, or export CLEO_PATH2CLEO.
###
### Steps (comma-separated, always done in this order; default: all):
###   configure  cmake configure
###   compile    compile the executable(s)
###   inputs     generate the input files
###   run        run the executable(s)
###   plot       plot the results
###   e.g. configure,compile or compile,run or plot
###
### To build an example from scratch, delete its build folder
### (by default <CLEO_PATH2CLEO>/build_xxx) before the configure step.
###
### Without an example, every entry in 'examples' below is used.
### An empty buildtype/compilername uses the machine default.
###
### Paths: edit the 'paths' section below before first use.
### Optional: debug build, make jobs and stack size in 'optional settings'.
### ============================================================ ###

set -e

### ------------- paths (EDIT THESE FOR YOUR SITE) ---------- ###
# Resolve paths from the submission directory while allowing site-specific overrides.
# Replace the <...> placeholders (or export the variables before running/submitting).
export CLEO_PATH2CLEO="${CLEO_PATH2CLEO:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
export CLEO_PYTHON="${CLEO_PYTHON:-${CLEO_PATH2CLEO}/.venv/bin/python3}"
export CLEO_YACYAXTROOT="${CLEO_YACYAXTROOT:-<PATH/TO/YACYAXT/INSTALL>}"
### -------------------------------------------------------- ###

### ----------- optional settings (defaults are fine) ------- ###
# Edit here, or export before running/submitting. Empty uses the machine default.
# exact build directory, only for a single example; empty builds each example
# in its default directory, <CLEO_PATH2CLEO>/build_xxx
export CLEO_PATH2BUILD="${CLEO_PATH2BUILD:-}"
export CLEO_ENABLEDEBUG="${CLEO_ENABLEDEBUG:-false}"   # true: debug build (-O0 -g)
export CLEO_MAKE_JOBS="${CLEO_MAKE_JOBS:-}"            # parallel make jobs
export CLEO_STACKSIZE_LIMIT="${CLEO_STACKSIZE_LIMIT:-}" # ulimit -s for runs (kB)
### -------------------------------------------------------- ###

### ---------------------- environment --------------------- ###
source /etc/profile
### -------------------------------------------------------- ###

### --------------------- configuration -------------------- ###
export CLEO_MACHINE="jupiter"

# "example buildtype compilername"
examples=(
  "as2017 cuda gcc"
  "cuspbifurc cuda gcc"
  "breakup cuda gcc"
  "shima2009 cuda gcc"
  "constthermo2d cuda gcc"
  "divfree2d cuda gcc"
  "eurec4a1d cuda gcc"
  "rainshaft1d cuda gcc"
  "python_bindings cuda gcc"
)
### -------------------------------------------------------- ###

source "${CLEO_PATH2CLEO}/scripts/common/run_jobs.sh"
# always a cuda build, whatever build type is given
run_cleo_jobs "${1:-}" "${2:-}" cuda "${4:-}"
