#!/bin/bash

### Core make logic for compiling CLEO executables.
### Requires CLEO_* environment variables to already be exported.
### Machine-specific module loading must be done by the calling wrapper BEFORE
### sourcing this script.

set -e

compile_cleo() {

  local executables=$1   # space-separated list (empty builds everything)
  local make_jobs=${CLEO_MAKE_JOBS:-8}  # override by setting CLEO_MAKE_JOBS (e.g. 128 on Levante)

  ### ---------------------------------------------------- ###

  ### ---------------- compile executables --------------- ###

  if [[ ! -f "${CLEO_PATH2BUILD}/Makefile" ]]; then
    echo "Error: no configured build found in ${CLEO_PATH2BUILD}."
    echo "Run the 'configure' step first."
    exit 1
  fi

  cd "${CLEO_PATH2BUILD}" && pwd

  echo "make -j ${make_jobs} ${executables}"
  make -j "${make_jobs}" ${executables}
  ### ---------------------------------------------------- ###
}

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  compile_cleo "$@"
fi
