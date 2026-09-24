#!/bin/bash

check_args_not_empty() {
  local inputs=("$@")
  for input in "${inputs[@]}"; do
    if [[ -z "$input" ]]; then
      echo "Bad inputs: please check all the required inputs have been specified"
      exit 1
    fi
  done
}

check_source_and_build_paths() {
  if [ "${CLEO_PATH2CLEO}" == "${CLEO_PATH2BUILD}" ]; then
    echo "Bad inputs: build directory cannot match the path to CLEO source"
    exit 1
  fi
}

check_buildtype() {
  if [[ "${CLEO_BUILDTYPE}" != "serial" &&
        "${CLEO_BUILDTYPE}" != "openmp" &&
        "${CLEO_BUILDTYPE}" != "threads" &&
        "${CLEO_BUILDTYPE}" != "cuda" ]];
  then
    echo "Bad inputs: build type must be 'serial', 'openmp', 'threads' or 'cuda'"
    exit 1
  fi
}

check_compilername() {
  if [[ "${CLEO_COMPILERNAME}" != "intel" && "${CLEO_COMPILERNAME}" != "gcc" ]]; then
    echo "Bad inputs: CLEO compiler name must be 'intel' or 'gcc'"
    exit 1
  fi
}

check_yac() {
  if [[ ${CLEO_YACYAXTROOT} == "" ]]
  then
    echo "Bad inputs: yacyaxtroot directory must be specified if YAC is enabled"
    exit 1
  fi
}
