#!/bin/bash
#SBATCH --job-name=build_compile_cleo
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --gpus=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=940M
#SBATCH --time=00:15:00
#SBATCH --mail-user=clara.bayley@mpimet.mpg.de
#SBATCH --mail-type=FAIL
#SBATCH --account=bm1183
#SBATCH --output=./build/bin/build_compile_cleo_out.%j.out
#SBATCH --error=./build/bin/build_compile_cleo_err.%j.out

set -e
source /etc/profile
module purge
spack unload --all

### ------------------ input parameters ---------------- ###
### ----- You need to edit these lines to specify ------ ###
### ----- your build configuration and executables ----- ###
### ---------------------------------------------------- ###
buildtype=$1                                   # "serial", "threads", "openmp" or "cuda"
compilername=${2:-intel}                       # "intel" or "gcc"
path2CLEO=${3:-${HOME}/CLEO}                   # must be absolute path
path2build=${4:-${path2CLEO}/build}            # should be absolute path
yacyaxtroot=${5:-/work/bm1183/m300950/yacyaxt/${compilername}} # yac and yaxt in yacyaxtroot/yac and yacyaxtroot/yaxt
build_flags=${6:-"-DCLEO_COUPLED_DYNAMICS="" \
  -DCLEO_PYTHON=/home/m/m300950/CLEO/.venv/bin/python3"} # CLEO_BUILD_FLAGS
executables=${7:-"cleocoupledsdm"}             # list of executables to compile or "NONE"
enabledebug=${8:-false}                        # == "true" or otherwise false
make_clean=${9:-true}                          # == "true" or otherwise false
### ---------------------------------------------------- ###

### ------------------ check arguments ----------------- ###
if [ "${path2CLEO}" == "" ]
then
  echo "Please provide path to CLEO source directory"
  exit 1
fi
source ${path2CLEO}/scripts/levante/bash/src/check_inputs.sh
check_args_not_empty "${buildtype}" "${path2CLEO}" "${path2build}"
check_args_not_empty "${compilername}" "${yacyaxtroot}" "${enabledebug}"
### ---------------------------------------------------- ###

### ----------------- export inputs -------------------- ###
export CLEO_BUILDTYPE=${buildtype}
export CLEO_COMPILERNAME=${compilername}
export CLEO_PATH2CLEO=${path2CLEO}
export CLEO_PATH2BUILD=${path2build}
export CLEO_BUILD_FLAGS=${build_flags}
export CLEO_YACYAXTROOT=${yacyaxtroot}
export CLEO_ENABLEDEBUG=${enabledebug}
### ---------------------------------------------------- ###

### -------------------- check inputs ------------------ ###
check_source_and_build_paths
check_buildtype
check_compilername
check_yac
### ---------------------------------------------------- ###

### -------------------- print inputs ------------------- ###
echo "### --------------- User Inputs -------------- ###"
echo "CLEO_BUILDTYPE = ${CLEO_BUILDTYPE}"
echo "CLEO_COMPILERNAME = ${CLEO_COMPILERNAME}"
echo "CLEO_PATH2CLEO = ${CLEO_PATH2CLEO}"
echo "CLEO_PATH2BUILD = ${CLEO_PATH2BUILD}"
echo "CLEO_BUILD_FLAGS = ${CLEO_BUILD_FLAGS}"
echo "CLEO_YACYAXTROOT = ${CLEO_YACYAXTROOT}"
echo "CLEO_ENABLEDEBUG = ${CLEO_ENABLEDEBUG}"
echo "executables = ${executables}"
echo "### ------------------------------------------- ###"
### ---------------------------------------------------- ###

### --------------------- build CLEO ------------------- ###
buildcmd="${CLEO_PATH2CLEO}/scripts/levante/bash/build_cleo.sh"
echo ${buildcmd}
eval ${buildcmd}
### ---------------------------------------------------- ###

### ---------------- compile executables --------------- ###
compilecmd="${CLEO_PATH2CLEO}/scripts/levante/bash/compile_cleo.sh \"${executables}\" ${make_clean}"
echo ${compilecmd}
eval ${compilecmd}
### ---------------------------------------------------- ###
