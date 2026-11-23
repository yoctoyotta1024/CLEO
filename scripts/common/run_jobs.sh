#!/bin/bash

### ============================================================ ###
###          Shared logic for the machine job scripts            ###
###          (<machine>/cpu.sh, <machine>/gpu.sh)                ###
### ============================================================ ###
###
### Usage (from a job script, after setting the configuration):
###   run_cleo_jobs [steps] [example] [buildtype] [compilername]
###
### Steps (comma-separated, always done in this order; default: all):
###   configure  cmake configure
###   compile    compile the executable(s)
###   inputs     generate the input files
###   run        run the executable(s)
###   plot       plot the results
###   e.g. configure,compile or compile,run or plot
###
### To build an example from scratch, delete its build folder
### (by default <CLEO_PATH2CLEO>/build_xxx) before the configure step.
###
### If an example is given only that one is used, otherwise every entry of
### the job script's 'examples' list ("example buildtype compilername").
###
### Job script configuration:
###   CLEO_MACHINE       machine name
###   CLEO_PATH2CLEO     CLEO source directory
###   CLEO_PYTHON        python used to run the examples
###   CLEO_YACYAXTROOT   YAC + YAXT installation
###   CLEO_PATH2BUILD    exact build directory, only for a single example (optional,
###                      default: each example in <CLEO_PATH2CLEO>/build_xxx)
###   CLEO_ENABLEDEBUG   true | false: debug build (optional, default: false)
###   CLEO_MAKE_JOBS     parallel make jobs (optional, default: machine)
###   CLEO_STACKSIZE_LIMIT  ulimit -s for runs in kB (optional, default: machine)
###   examples=(...)  default list of "example buildtype compilername"
### ============================================================ ###

set -e

run_cleo_jobs() {
  local steps="${1:-all}"
  local example="${2:-}"
  local buildtype="${3:-}"
  local compilername="${4:-}"

  source "${CLEO_PATH2CLEO}/scripts/common/check_inputs.sh"
  check_steps "${steps}"
  check_machine

  check_args_not_empty "${CLEO_PATH2CLEO}" "${CLEO_PYTHON}" "${CLEO_YACYAXTROOT}"

  # catch paths that were left as <...> placeholders in the job script
  local var
  for var in CLEO_PATH2CLEO CLEO_PYTHON CLEO_YACYAXTROOT CLEO_PATH2BUILD; do
    if [[ "${!var}" == *"<"*">"* ]]; then
      echo "Error: ${var} is still a placeholder: ${!var}"
      echo "Edit the 'paths' section of $0 or export ${var} before running."
      exit 1
    fi
  done
  if [[ "${steps}" == all || ",${steps}," =~ ,(inputs|run|plot), ]] &&
     [[ ! -x "${CLEO_PYTHON}" ]] && ! command -v "${CLEO_PYTHON}" &>/dev/null; then
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

  # one exact build directory would be shared (and overwritten) by every example
  if [[ -n "${CLEO_PATH2BUILD:-}" && ${#entries[@]} -gt 1 ]]; then
    echo "Error: CLEO_PATH2BUILD is set, so give a single example to run."
    echo "Leave CLEO_PATH2BUILD empty to build each example in its default directory."
    exit 1
  fi

  local entry e b c
  for entry in "${entries[@]}"; do
    read -r e b c <<< "${entry}"
    [[ "${b}" == "-" ]] && b=""
    [[ "${c}" == "-" ]] && c=""
    echo
    echo "=== ${steps}: ${e} (${b:-default buildtype}, ${c:-default compiler}) ==="
    echo
    "${driver}" "${e}" "${b}" "${c}" "${CLEO_PATH2CLEO}" \
      "${CLEO_PATH2BUILD}" "" "${CLEO_YACYAXTROOT}" "" "" "${steps}"
  done
}
