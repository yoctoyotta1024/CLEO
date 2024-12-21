#!/bin/bash
#SBATCH --job-name=compile_cleo
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=940M
#SBATCH --time=00:05:00
#SBATCH --mail-user=clara.bayley@mpimet.mpg.de
#SBATCH --mail-type=FAIL
#SBATCH --account=bm1183
#SBATCH --output=./build/bin/compile_cleo_out.%j.out
#SBATCH --error=./build/bin/compile_cleo_err.%j.out

### Please note: script may assume required CLEO_[XXX]
### variables have already exported (!)

set -e
module purge
spack unload --all

executables=$1
make_clean=$2
bashsrc=${CLEO_PATH2CLEO}/scripts/bash/src

### -------------------- check inputs ------------------ ###
source ${bashsrc}/check_inputs.sh
check_args_not_empty "${CLEO_BUILDTYPE}" "${CLEO_COMPILERNAME}" \
  "${CLEO_PATH2CLEO}" "${CLEO_PATH2BUILD}"

check_source_and_build_paths
check_buildtype
check_compilername
### ---------------------------------------------------- ###

### ----------------- load compiler(s) ----------------- ###
if [ "${CLEO_COMPILERNAME}" == "intel" ]
then
  echo "TODO(CB): intel compiler support"
  exit 1
elif [ "${CLEO_COMPILERNAME}" == "gcc" ]
then
  module load gcc/11.2.0-gcc-11.2.0 openmpi/4.1.2-gcc-11.2.0
  spack load cmake@3.23.1%gcc
  if [ "${CLEO_BUILDTYPE}" == "cuda" ]
  then
    module load nvhpc/23.9-gcc-11.2.0
  fi
  echo "TODO(CB): update gcc compiler version (in YAC and cuda too!)"
  exit 1
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
cmd="make -j 128 ${executables}"
echo ${cmd}
eval ${cmd}
### ---------------------------------------------------- ###
