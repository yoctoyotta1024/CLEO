#!/bin/bash

### ============================================================ ###
###                    Vanilla CPU job script                    ###
### ============================================================ ###
###
### Usage:
###   ./scripts/vanilla/cpu.sh build [experiment] [buildtype] [compilername]
###   ./scripts/vanilla/cpu.sh [all|run] [experiment] [buildtype] [compilername]
###
###   Run from the CLEO root directory.
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
# no module system: mpic++, mpicc and cmake must already be on PATH
### -------------------------------------------------------- ###

### --------------------- configuration -------------------- ###
export CLEO_MACHINE="vanilla"

# "experiment buildtype compilername"
experiments=(
  "as2017 serial gcc"
)

# true: delete each experiment's build folder first and rebuild from scratch
export CLEO_MAKE_CLEAN="${CLEO_MAKE_CLEAN:-false}"

# command prefix for run mode (e.g. srun), empty to run directly
run_launcher=()
### -------------------------------------------------------- ###

source "${CLEO_PATH2CLEO}/scripts/common/run_jobs.sh"
run_cleo_jobs "$@"
