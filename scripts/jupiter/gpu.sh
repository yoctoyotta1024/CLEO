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

### ============================================================ ###
###                    Jupiter GPU job script                    ###
### ============================================================ ###
###
### Usage:
###   ./scripts/jupiter/gpu.sh build [experiment] [buildtype] [compilername]
###   sbatch scripts/jupiter/gpu.sh [all|run] [experiment] [buildtype] [compilername]
###
###   Run from (or submit from) the CLEO root directory.
###
### Modes:
###   all    configure, compile, run + plot (default)    (steps: all)
###   build  cmake configure + compile                   (steps: build,compile)
###   run    recompile, run + plot (needs a prior build) (steps: compile,run,plot)
###
### Set CLEO_MAKE_CLEAN=true below (or export it) to delete each
### experiment's build folder and rebuild from scratch (all and build modes).
###
### Without an experiment, every entry in 'experiments' below is used.
### An empty buildtype/compilername uses the machine default.
###
### Paths: edit the 'paths' section below before first use.
### ============================================================ ###

set -e

### ------------- paths (EDIT THESE FOR YOUR SITE) ---------- ###
# Resolve paths from the submission directory while allowing site-specific overrides.
# Replace the <...> placeholders (or export the variables before running/submitting).
# CLEO_PATH2BUILD is a build root: each experiment builds in <CLEO_PATH2BUILD>/build_xxx
export CLEO_PATH2CLEO="${SLURM_SUBMIT_DIR:-$(pwd)}"
export CLEO_PYTHON="${CLEO_PYTHON:-${CLEO_PATH2CLEO}/.venv/bin/python3}"
export CLEO_YACYAXTROOT="${CLEO_YACYAXTROOT:-<PATH/TO/YACYAXT/INSTALL>}"
export CLEO_PATH2BUILD="${CLEO_PATH2BUILD:-<PATH/TO/BUILD/ROOT>}"
### -------------------------------------------------------- ###

### ---------------------- environment --------------------- ###
source /etc/profile
# Required if the .venv does not use the default environment.
module load Stages/2026 Python/3.13.5
### -------------------------------------------------------- ###

### --------------------- configuration -------------------- ###
export CLEO_MACHINE="jupiter"

# "experiment buildtype compilername"
experiments=(
  "constthermo2d cuda gcc"
)

# true: delete each experiment's build folder first and rebuild from scratch
export CLEO_MAKE_CLEAN="${CLEO_MAKE_CLEAN:-false}"
### -------------------------------------------------------- ###

source "${CLEO_PATH2CLEO}/scripts/common/run_jobs.sh"
run_cleo_jobs "$@"
