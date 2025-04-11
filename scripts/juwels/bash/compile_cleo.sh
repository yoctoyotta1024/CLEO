#!/bin/bash
#SBATCH --job-name=compile_cleo
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=940M
#SBATCH --time=00:05:00
#SBATCH --account=exaww
#SBATCH --output=./build/bin/compile_cleo_out.%j.out
#SBATCH --error=./build/bin/compile_cleo_err.%j.out

### Please note: script may assume required CLEO_[XXX]
### variables have already exported (!)

set -e
module purge

executables=$1                     # if == "NONE" only libraries built
make_clean=$2

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
bashsrc=${SCRIPT_DIR}/src

### -------------------- check inputs ------------------ ###
source ${bashsrc}/check_inputs.sh
check_args_not_empty "${CLEO_BUILDTYPE}" "${CLEO_COMPILERNAME}" \
  "${CLEO_PATH2CLEO}" "${CLEO_PATH2BUILD}"

check_source_and_build_paths
check_buildtype
check_compilername
### ---------------------------------------------------- ###

### ----------------- load compiler(s) ----------------- ###
source ${bashsrc}/juwels_packages.sh

if [ "${CLEO_COMPILERNAME}" == "intel" ]
then
  module load ${juwels_intel}
  module load ${juwels_intel_mpi}
  module load ${juwels_intel_cmake}
elif [ "${CLEO_COMPILERNAME}" == "gcc" ]
then
  module load ${juwels_gcc}
  module load ${juwels_gcc_mpi}
  module load ${juwels_gcc_cmake}
  if [ "${CLEO_BUILDTYPE}" == "cuda" ]
  then
    echo "Bad inputs, CUDA build enabled but building CLEO with CUDA on JUWELS is not currently supported"
    exit 1
  fi
fi
### ---------------------------------------------------- ###

### ---------------- compile executables --------------- ###
echo "### --------------- Compile Inputs -------------- ###"
echo "CLEO_BUILDTYPE: ${CLEO_BUILDTYPE}"
echo "CLEO_COMPILERNAME: ${CLEO_COMPILERNAME}"
echo "CLEO_PATH2CLEO: ${CLEO_PATH2CLEO}"
echo "CLEO_PATH2BUILD: ${CLEO_PATH2BUILD}"

echo "executables: ${executables}"
echo "make_clean: ${make_clean}"
echo "### ------------------------------------------- ###"

cd ${CLEO_PATH2BUILD} && pwd

if [ "${make_clean}" == "true" ]
then
  cmd="make clean"
  echo ${cmd}
  eval ${cmd}
fi

if [ ${executables} == "NONE" ]
then
  cmd="make -j 96"
else
  cmd="make -j 96 ${executables}"
fi
echo ${cmd}
eval ${cmd}
### ---------------------------------------------------- ###
