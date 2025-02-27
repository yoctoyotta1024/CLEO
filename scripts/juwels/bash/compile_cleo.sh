#!/bin/bash
#SBATCH --job-name=compile_cleo
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=940M
#SBATCH --time=00:05:00
#SBATCH --mail-user=clara.bayley@mpimet.mpg.de
#SBATCH --mail-type=FAIL
#SBATCH --account=exaww
#SBATCH --output=./build/bin/compile_cleo_out.%j.out
#SBATCH --error=./build/bin/compile_cleo_err.%j.out

### Please note: script may assume required CLEO_[XXX]
### variables have already exported (!)

set -e
module purge

executables=$1
make_clean=$2
bashsrc=${CLEO_PATH2CLEO}/scripts/juwels/bash/src

### -------------------- check inputs ------------------ ###
source ${bashsrc}/check_inputs.sh
check_args_not_empty "${CLEO_BUILDTYPE}" "${CLEO_COMPILERNAME}" \
  "${CLEO_PATH2CLEO}" "${CLEO_PATH2BUILD}" "${CLEO_ENABLE_MPTRAC}"

check_source_and_build_paths
check_buildtype
check_compilername
### ---------------------------------------------------- ###


### ----------------- load compiler(s) ----------------- ###
bashsrc=${CLEO_PATH2CLEO}/scripts/juwels/bash/src
source ${bashsrc}/juwels_packages.sh

if [ "${CLEO_BUILDTYPE}" == "cuda" ]
then
  echo "Bad inputs, CUDA build enabled but building CLEO with CUDA on JUWELS is not currently supported"
  exit 1
fi

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
fi
### ---------------------------------------------------- ###

### ----------- MPTRAC compile-time settings ----------- ###
if [ "${CLEO_ENABLE_MPTRAC}" == "true" ]
then
  if [ "${CLEO_COMPILERNAME}" != "gcc" ]
  then
    echo "Bad inputs, MPTRAC can only supported with gcc compilation on JUWELS"
    exit 1
  fi
  module load ${juwels_gcc_gsl}
  module load ${juwels_gcc_netcdf}
fi
### ---------------------------------------------------- ###

### ---------------- compile executables --------------- ###
echo "### --------------- Compile Inputs -------------- ###"
echo "CLEO_BUILDTYPE: ${CLEO_BUILDTYPE}"
echo "CLEO_COMPILERNAME: ${CLEO_COMPILERNAME}"
echo "CLEO_PATH2CLEO: ${CLEO_PATH2CLEO}"
echo "CLEO_PATH2BUILD: ${CLEO_PATH2BUILD}"
echo "CLEO_ENABLE_MPTRAC: ${CLEO_ENABLE_MPTRAC}"

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
cmd="make -j 96 ${executables}"
echo ${cmd}
eval ${cmd}
### ---------------------------------------------------- ###
