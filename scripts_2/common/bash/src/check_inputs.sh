#!/bin/bash

### Shared input validation functions used across all machines.
### Valid build types / compiler names are the union of all supported machines.

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

check_yac() {
  if [[ ${CLEO_YACYAXTROOT} == "" ]]; then
    echo "Bad inputs: yacyaxtroot directory must be specified if YAC is enabled"
    exit 1
  fi
}

check_machine() {
  if [[ -z "${CLEO_MACHINE}" ]]; then
    echo "Bad inputs: CLEO_MACHINE must be set (e.g. 'vanilla', 'levante', 'juwels')"
    exit 1
  fi
}
