#!/bin/bash

### ============================================================ ###
###                        Jupiter script                        ###
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
###   buildtype        serial threads openmp cuda (default: openmp)
###   compilername     gcc (default: gcc)
###   stacksize_limit  204800 (kB)
###   CLEO_MAKE_JOBS   32
###
###   steps            build,compile,run,plot,all (default: all)
###                    The run stage generates input files and runs the executable.
###                    Set NSYS_PREFIX to profile only the executable run.
###
### Supported experiments (see common/experiments.sh for full details):
###   as2017 cuspbifurc breakup shima2009 constthermo2d divfree2d
###   eurec4a1d rainshaft1d python_bindings fromfile
###   fromfile_irreg bubble3d
### ============================================================ ###

set -e

### ---------------- machine configuration ------------- ###
export CLEO_MACHINE="jupiter"
machine_default_buildtype="openmp"
machine_buildtypes=(serial threads openmp cuda)
machine_compilers=(gcc)
machine_experiments=(
  as2017 cuspbifurc breakup shima2009 constthermo2d divfree2d
  eurec4a1d rainshaft1d python_bindings fromfile
  fromfile_irreg bubble3d
)
machine_default_stacksize="204800"
machine_default_make_jobs=32
### ---------------------------------------------------- ###

source "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../common" &>/dev/null && pwd)/build_compile_run_plot_cleo.sh"
build_compile_run_plot_cleo "$@"
