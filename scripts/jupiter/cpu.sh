#!/bin/bash
#SBATCH --job-name=cleo_cpu
#SBATCH --partition=booster
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=288
#SBATCH --time=00:30:00
#SBATCH --account=xspies
#SBATCH --output=./cleo_cpu.%j.out
#SBATCH --error=./cleo_cpu.%j.out

### ============================================================ ###
###                    Jupiter CPU job script                    ###
### ============================================================ ###
###
### Usage:
###   ./scripts/jupiter/cpu.sh build [example] [buildtype] [compilername]
###   sbatch scripts/jupiter/cpu.sh [all|run] [example] [buildtype] [compilername]
###
###   Run from (or submit from) the CLEO root directory.
###
### Modes:
###   all    configure, compile, run + plot (default)    (steps: all)
###   build  cmake configure + compile                   (steps: build,compile)
###   run    recompile, run + plot (needs a prior build) (steps: compile,run,plot)
###
### Set CLEO_MAKE_CLEAN=true below (or export it) to delete each
### example's build folder and rebuild from scratch (all and build modes).
###
### Without an example, every entry in 'examples' below is used.
### An empty buildtype/compilername uses the machine default.
###
### Paths: edit the 'paths' section below before first use.
### ============================================================ ###

set -e

### ------------- paths (EDIT THESE FOR YOUR SITE) ---------- ###
# Resolve paths from the submission directory while allowing site-specific overrides.
# Replace the <...> placeholders (or export the variables before running/submitting).
# CLEO_PATH2BUILD is a build root: each example builds in <CLEO_PATH2BUILD>/build_xxx
export CLEO_PATH2CLEO="${SLURM_SUBMIT_DIR:-$(pwd)}"
export CLEO_PYTHON="${CLEO_PYTHON:-${CLEO_PATH2CLEO}/.venv/bin/python3}"
export CLEO_YACYAXTROOT="${CLEO_YACYAXTROOT:-<PATH/TO/YACYAXT/INSTALL>}"
export CLEO_PATH2BUILD="${CLEO_PATH2BUILD:-<PATH/TO/BUILD/ROOT>}"
### -------------------------------------------------------- ###

### ---------------------- environment --------------------- ###
source /etc/profile
### -------------------------------------------------------- ###

### --------------------- configuration -------------------- ###
export CLEO_MACHINE="jupiter"

# "example buildtype compilername"
examples=(
  "constthermo2d openmp gcc"
)

# true: delete each example's build folder first and rebuild from scratch
export CLEO_MAKE_CLEAN="${CLEO_MAKE_CLEAN:-false}"
### -------------------------------------------------------- ###

source "${CLEO_PATH2CLEO}/scripts/common/run_jobs.sh"
run_cleo_jobs "$@"
