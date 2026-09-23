#!/bin/bash

### ============================================================ ###
###                        Vanilla script                        ###
### ============================================================ ###
###
### Usage:
###   ./build_compile_run_plot_cleo.sh [experiment] [buildtype] [compilername] \
###                                    [path2CLEO] [path2build] [build_flags] \
###                                    [yacyaxtroot] [enabledebug] [make_clean] \
###                                    [stacksize_limit] [steps]
###
###   All arguments are optional; pass "" to keep a default.
###   path2build is a build root: the experiment builds in <path2build>/build_xxx.
###   The arguments are the same on every machine, see
###   common/build_compile_run_plot_cleo.sh for their full description.
###
### Machine defaults / supported values:
###   buildtype        serial threads openmp (default: serial)
###   compilername     gcc (default: gcc)
###   stacksize_limit  unchanged (kB)
###   CLEO_MAKE_JOBS   8
###
###   make_clean       true deletes the build folder first (build from scratch).
###   steps            build,compile,run,plot,all (default: all)
###                    The run stage generates input files and runs the executable.
###                    Set NSYS_PREFIX to profile only the executable run.
###
### Supported experiments (see common/experiments.sh for full details):
###   as2017 cuspbifurc breakup shima2009 constthermo2d divfree2d
###   eurec4a1d rainshaft1d python_bindings
###   Note: fromfile, fromfile_irreg and bubble3d need an HPC machine
###   (levante or jupiter).
### ============================================================ ###

set -e

### ---------------- machine configuration ------------- ###
export CLEO_MACHINE="vanilla"
machine_default_buildtype="serial"
machine_buildtypes=(serial threads openmp)
machine_compilers=(gcc)
machine_experiments=(
  as2017 cuspbifurc breakup shima2009 constthermo2d divfree2d
  eurec4a1d rainshaft1d python_bindings
)
machine_default_stacksize=""
machine_default_make_jobs=8
### ---------------------------------------------------- ###

source "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../common" &>/dev/null && pwd)/build_compile_run_plot_cleo.sh"
build_compile_run_plot_cleo "$@"
