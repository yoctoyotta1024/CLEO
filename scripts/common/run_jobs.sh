#!/bin/bash

### ============================================================ ###
###          Shared logic for the machine job scripts            ###
###          (<machine>/cpu.sh, <machine>/gpu.sh)                ###
### ============================================================ ###
###
### Usage (from a job script, after setting the configuration):
###   run_cleo_jobs [all|build|run] [example] [buildtype] [compilername]
###
### Modes:
###   all    configure, compile, run + plot (default)    (steps: all)
###   build  cmake configure + compile                   (steps: build,compile)
###   run    recompile, make inputs, run + plot (needs a prior build)
###          (steps: compile,inputs,run,plot)
###
### To build an example from scratch, delete its build folder
### (<CLEO_PATH2BUILD>/build_xxx) before using the all or build mode.
###
### If an example is given only that one is used, otherwise every entry of
### the job script's 'examples' list ("example buildtype compilername").
###
### Job script configuration:
###   CLEO_MACHINE       machine name
###   CLEO_PATH2CLEO     CLEO source directory
###   CLEO_PYTHON        python used to run the examples
###   CLEO_YACYAXTROOT   YAC + YAXT installation
###   CLEO_PATH2BUILD    build root: examples build in <CLEO_PATH2BUILD>/build_xxx
###   CLEO_ENABLEDEBUG   true | false: debug build (optional, default: false)
###   CLEO_MAKE_JOBS     parallel make jobs (optional, default: machine)
###   CLEO_STACKSIZE_LIMIT  ulimit -s for runs in kB (optional, default: machine)
###   examples=(...)  default list of "example buildtype compilername"
### ============================================================ ###

set -e

run_cleo_jobs() {
  local mode="${1:-all}"
  local example="${2:-}"
  local buildtype="${3:-}"
  local compilername="${4:-}"

  local steps
  local action
  case "${mode}" in
    build)
      steps="build,compile"
      action="Building"
      ;;
    run)
      steps="compile,inputs,run,plot"
      action="Running"
      ;;
    all)
      steps="all"
      action="Building + running"
      ;;
    *)
      echo "Usage: $0 [all|build|run] [example] [buildtype] [compilername]"
      exit 1
      ;;
  esac

  source "${CLEO_PATH2CLEO}/scripts/common/check_inputs.sh"
  check_machine

  check_args_not_empty "${CLEO_PATH2CLEO}" "${CLEO_PYTHON}" "${CLEO_YACYAXTROOT}" "${CLEO_PATH2BUILD}"

  # catch paths that were left as <...> placeholders in the job script
  local var
  for var in CLEO_PATH2CLEO CLEO_PYTHON CLEO_YACYAXTROOT CLEO_PATH2BUILD; do
    if [[ "${!var}" == *"<"*">"* ]]; then
      echo "Error: ${var} is still a placeholder: ${!var}"
      echo "Edit the 'paths' section of $0 or export ${var} before running."
      exit 1
    fi
  done
  if [[ "${mode}" != build && ! -x "${CLEO_PYTHON}" ]] && ! command -v "${CLEO_PYTHON}" &>/dev/null; then
    echo "Error: CLEO Python executable not found:"
    echo "  ${CLEO_PYTHON}"
    exit 1
  fi

  local driver="${CLEO_PATH2CLEO}/scripts/${CLEO_MACHINE}/build_compile_run_plot_cleo.sh"

  local entries=()
  if [[ -n "${example}" ]]; then
    entries=("${example} ${buildtype:--} ${compilername:--}")
  else
    entries=("${examples[@]}")
  fi

  local entry e b c
  for entry in "${entries[@]}"; do
    read -r e b c <<< "${entry}"
    [[ "${b}" == "-" ]] && b=""
    [[ "${c}" == "-" ]] && c=""
    echo
    echo "=== ${action} ${e} (${b:-default buildtype}, ${c:-default compiler}) ==="
    echo
    "${driver}" "${e}" "${b}" "${c}" "${CLEO_PATH2CLEO}" \
      "${CLEO_PATH2BUILD}" "" "${CLEO_YACYAXTROOT}" "" "" "${steps}"
  done
}
