#!/bin/bash

### Shared input validation functions used across all machines.

check_args_not_empty() {
  local inputs=("$@")
  for input in "${inputs[@]}"; do
    if [[ -z "$input" ]]; then
      echo "Bad inputs: please check all the required inputs have been specified"
      exit 1
    fi
  done
}

# usage: check_value_in_list <name> <value> <allowed values...>
check_value_in_list() {
  local name="$1"
  local value="$2"
  shift 2
  local allowed
  for allowed in "$@"; do
    [[ "${value}" == "${allowed}" ]] && return 0
  done
  echo "Bad inputs: ${name} '${value}' is not supported on ${CLEO_MACHINE:-this machine}."
  echo "Supported: $*"
  exit 1
}

# usage: check_steps <steps>  (comma-separated build,compile,run,plot, or all)
check_steps() {
  local steps="$1"
  local step
  local requested_steps
  [[ "${steps}" == all ]] && return 0
  IFS=',' read -r -a requested_steps <<< "${steps}"
  for step in "${requested_steps[@]}"; do
    case "${step}" in
      build|compile|run|plot) ;;
      *)
        echo "Bad inputs: invalid step '${step}'. Use build, compile, run, plot, or all."
        exit 1
        ;;
    esac
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
    echo "Bad inputs: CLEO_MACHINE must be set (e.g. 'vanilla', 'levante', 'jupiter')"
    exit 1
  fi
}
