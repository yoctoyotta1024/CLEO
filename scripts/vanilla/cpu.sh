#!/bin/bash

### ============================================================ ###
###                    Vanilla CPU job script                    ###
### ============================================================ ###
###
### Usage:
###   ./scripts/vanilla/cpu.sh build [example] [buildtype] [compilername]
###   ./scripts/vanilla/cpu.sh [all|run] [example] [buildtype] [compilername]
###
###   Run from the CLEO root directory.
###
### Modes:
###   all    configure, compile, run + plot (default)    (steps: all)
###   build  cmake configure + compile                   (steps: build,compile)
###   run    recompile, run + plot (needs a prior build) (steps: compile,run,plot)
###
### To build an example from scratch, delete its build folder
### (<CLEO_PATH2BUILD>/build_xxx) before using the all or build mode.
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
# CLEO_PATH2BUILD is a build root: each example builds in <CLEO_PATH2BUILD>/build_xxx
export CLEO_PATH2CLEO="${SLURM_SUBMIT_DIR:-$(pwd)}"
export CLEO_PYTHON="${CLEO_PYTHON:-${CLEO_PATH2CLEO}/.venv/bin/python3}"
export CLEO_YACYAXTROOT="${CLEO_YACYAXTROOT:-<PATH/TO/YACYAXT/INSTALL>}"
export CLEO_PATH2BUILD="${CLEO_PATH2BUILD:-<PATH/TO/BUILD/ROOT>}"
### -------------------------------------------------------- ###

### ----------- optional settings (defaults are fine) ------- ###
# Edit here, or export before running/submitting. Empty uses the machine default.
export CLEO_ENABLEDEBUG="${CLEO_ENABLEDEBUG:-false}"   # true: debug build (-O0 -g)
export CLEO_MAKE_JOBS="${CLEO_MAKE_JOBS:-}"            # parallel make jobs
export CLEO_STACKSIZE_LIMIT="${CLEO_STACKSIZE_LIMIT:-}" # ulimit -s for runs (kB)
### -------------------------------------------------------- ###

### ---------------------- environment --------------------- ###
# no module system: mpic++, mpicc and cmake must already be on PATH
### -------------------------------------------------------- ###

### --------------------- configuration -------------------- ###
export CLEO_MACHINE="vanilla"

# "example buildtype compilername"
examples=(
  "as2017 serial gcc"
  "cuspbifurc openmp gcc"
  "breakup serial gcc"
  "shima2009 openmp gcc"
  "constthermo2d serial gcc"
  "divfree2d openmp gcc"
  "eurec4a1d serial gcc"
  "rainshaft1d openmp gcc"
  # known to fail at the moment, so last (the job stops at the first failure):
  "roughpaper serial gcc"
  "python_bindings openmp gcc"
)
### -------------------------------------------------------- ###

source "${CLEO_PATH2CLEO}/scripts/common/run_jobs.sh"
run_cleo_jobs "$@"
