# Usage: print_configuration.sh <EXPERIMENT>
#
# Prints the current build configuration to stdout.
#
# Parameters:
#   $1  EXPERIMENT  - Name or identifier of the experiment to run.
#

print_configuration() {
  local experiment_name=$1

  echo "### --------------- User Inputs -------------- ###"
  echo "CLEO_BUILDTYPE = ${CLEO_BUILDTYPE}"
  echo "CLEO_COMPILERNAME = ${CLEO_COMPILERNAME}"
  echo "CLEO_PATH2CLEO = ${CLEO_PATH2CLEO}"
  echo "CLEO_PATH2BUILD = ${CLEO_PATH2BUILD}"
  echo "CLEO_BUILD_FLAGS = ${CLEO_BUILD_FLAGS}"
  echo "CLEO_YACYAXTROOT = ${CLEO_YACYAXTROOT}"
  echo "CLEO_ENABLEDEBUG = ${CLEO_ENABLEDEBUG}"
  echo "EXPERIMENT = ${experiment_name}"
  echo "### ------------------------------------------- ###"
}

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  print_configuration "$@"
fi
