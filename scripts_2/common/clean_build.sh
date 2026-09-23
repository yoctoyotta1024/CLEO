#!/bin/bash

### Deletes an experiment's build folder so the next build starts from scratch.
### Requires CLEO_PATH2BUILD and CLEO_PATH2CLEO to be exported.
###
### For safety it only deletes a folder that contains a CMakeCache.txt, and
### never /, $HOME, the CLEO source folder or a folder containing it.

set -e

clean_cleo_build() {
  local common_dir
  common_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)
  source "${common_dir}/check_inputs.sh"
  check_args_not_empty "${CLEO_PATH2BUILD}" "${CLEO_PATH2CLEO}"

  local build_dir="${CLEO_PATH2BUILD}"
  if [[ ! -e "${build_dir}" ]]; then
    echo "Nothing to clean: ${build_dir} does not exist"
    return 0
  fi

  local resolved cleo_src home
  resolved=$(cd -- "${build_dir}" && pwd -P)
  cleo_src=$(cd -- "${CLEO_PATH2CLEO}" && pwd -P)
  home=$(cd -- "${HOME}" && pwd -P)

  if [[ "${resolved}" == "/" || "${resolved}" == "${home}" || \
        "${resolved}" == "${cleo_src}" || "${cleo_src}" == "${resolved}/"* ]]; then
    echo "Error: refusing to delete ${resolved} (it is /, \$HOME or contains the CLEO source)."
    exit 1
  fi

  if [[ ! -f "${resolved}/CMakeCache.txt" ]]; then
    echo "Error: refusing to delete ${resolved}: no CMakeCache.txt, so it does not look like a CLEO build folder."
    echo "Delete it by hand if you are sure."
    exit 1
  fi

  echo "Removing build folder to build from scratch: ${resolved}"
  rm -rf -- "${resolved}"
}

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  clean_cleo_build "$@"
fi
