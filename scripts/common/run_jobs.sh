#!/bin/bash

### ============================================================ ###
###          Shared logic for the machine job scripts            ###
###          (<machine>/cpu.sh, <machine>/gpu.sh)                ###
### ============================================================ ###
###
### Usage (from a job script, after setting the configuration):
###   run_cleo_jobs [all|build|run] [experiment] [buildtype] [compilername]
###
### Modes:
###   all    configure, compile, run + plot (default)    (steps: all)
###   build  cmake configure + compile                   (steps: build,compile)
###   run    recompile, run + plot (needs a prior build) (steps: compile,run,plot)
###
### Set CLEO_MAKE_CLEAN=true to delete each experiment's build folder before
### building (build and all modes only).
###
### If an experiment is given only that one is used, otherwise every entry of
### the job script's 'experiments' list ("experiment buildtype compilername").
###
### Job script configuration:
###   CLEO_MACHINE       machine name
###   CLEO_PATH2CLEO     CLEO source directory
###   CLEO_PYTHON        python used to run the experiments
###   CLEO_YACYAXTROOT   YAC + YAXT installation
###   CLEO_PATH2BUILD    build root: experiments build in <CLEO_PATH2BUILD>/build_xxx
###   CLEO_MAKE_CLEAN    true | false: delete build folders first (default: false)
###   experiments=(...)  default list of "experiment buildtype compilername"
### ============================================================ ###

set -e

run_cleo_jobs() {
  local mode="${1:-all}"
  local experiment="${2:-}"
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
      steps="compile,run,plot"
      action="Running"
      ;;
    all)
      steps="all"
      action="Building + running"
      ;;
    *)
      echo "Usage: $0 [all|build|run] [experiment] [buildtype] [compilername]"
      exit 1
      ;;
  esac

  source "${CLEO_PATH2CLEO}/scripts/common/check_inputs.sh"
  check_machine

  local make_clean="${CLEO_MAKE_CLEAN:-false}"
  check_value_in_list CLEO_MAKE_CLEAN "${make_clean}" true false
  if [[ "${make_clean}" == true ]]; then
    if [[ "${mode}" == run ]]; then
      echo "Error: CLEO_MAKE_CLEAN=true needs 'all' or 'build' mode ('run' reuses an existing build)."
      exit 1
    fi
    action="${action} from scratch"
  fi
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
  if [[ -n "${experiment}" ]]; then
    entries=("${experiment} ${buildtype:--} ${compilername:--}")
  else
    entries=("${experiments[@]}")
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
      "${CLEO_PATH2BUILD}" "" "${CLEO_YACYAXTROOT}" "" "${make_clean}" "" "${steps}"
  done
}
