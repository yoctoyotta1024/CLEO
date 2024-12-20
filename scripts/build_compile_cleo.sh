#!/bin/bash
#SBATCH --job-name=build_compile_cleo
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=940M
#SBATCH --time=00:05:00
#SBATCH --mail-user=clara.bayley@mpimet.mpg.de
#SBATCH --mail-type=FAIL
#SBATCH --account=bm1183
#SBATCH --output=./build/bin/build_compile_cleo_out.%j.out
#SBATCH --error=./build/bin/build_compile_cleo_err.%j.out

set -e
module purge
spack unload --all

### ------------------ input parameters ---------------- ###
### ----- You need to edit these lines to specify ------ ###
### ----- your build configuration and executables ----- ###
### ---------------------------------------------------- ###
buildtype=$1                                # "serial", "threads", "openmp" or "cuda"
compilername=${2:-intel}                    # "intel" or "gcc"
path2CLEO=${3:-${HOME}/CLEO}                # must be absolute path
path2build=${4:-${path2CLEO}/build}         # should be absolute path
enabledebug=${5:-false}                     # "true" or otherwise false
enableyac=${6:-false}                       # "true" or otherwise false
executables=${7:-"cleocoupledsdm testing"}  # list of executables to compile
yacyaxtroot=/work/bm1183/m300950/yacyaxt    # yac and yaxt in yacyaxtroot/yac and yacyaxtroot/yaxt
### ---------------------------------------------------- ###

### -------------------- check inputs ------------------ ###
if [[ "${buildtype}" == "" || "${compilername}" == "" || "${enabledebug}" == "" ||
      "${path2CLEO}" == "" || "${path2build}" == "" || "${enableyac}" == "" ]]
then
  echo "Bad inputs, please check all the required inputs have been specified"
  exit 1
fi

if [[ "${path2CLEO}" == "${path2build}" ]]
then
  echo "Bad inputs, build directory cannot match the path to CLEO source"
  exit 1
fi

if [ "${CLEO_BUILDTYPE}" != "serial" ] &&
   [ "${CLEO_BUILDTYPE}" != "openmp" ] &&
   [ "${CLEO_BUILDTYPE}" != "threads" ] &&
   [ "${CLEO_BUILDTYPE}" != "cuda" ];
then
  echo "Bad inputs, build type must be 'serial', 'openmp', 'threads' or 'cuda'"
  exit 1
fi

if [ ${CLEO_ENABLEYAC} == "true" && ${CLEO_YACYAXTROOT} == "" ]
then
  echo "Bad inputs, yacyaxtroot directory must be specified if YAC is enabled"
  exit 1
fi
### ---------------------------------------------------- ###

### ----------------- export inputs -------------------- ###
export CLEO_BUILDTYPE=${buildtype}
export CLEO_COMPILERNAME=${compilername}
export CLEO_PATH2CLEO=${path2CLEO}
export CLEO_PATH2BUILD=${path2build}
export CLEO_ENABLEDEBUG=${enabledebug}
export CLEO_ENABLEYAC=${enableyac}

if [ ${CLEO_ENABLEYAC} == "true" ]
then
  export CLEO_YACYAXTROOT=${yacyaxtroot}
fi
### ---------------------------------------------------- ###

### -------------------- print inputs ------------------- ###
echo "### --------------- User Inputs -------------- ###"
echo "CLEO_BUILDTYPE = ${CLEO_BUILDTYPE}"
echo "CLEO_COMPILERNAME = ${CLEO_COMPILERTYPE}"
echo "CLEO_PATH2CLEO = ${CLEO_PATH2CLEO}"
echo "CLEO_PATH2BUILD = ${CLEO_PATH2BUILD}"
echo "CLEO_ENABLEDEBUG = ${CLEO_ENABLEDEBUG}"
echo "CLEO_ENABLEYAC = ${CLEO_ENABLEYAC}"
echo "CLEO_YACYAXTROOT = ${CLEO_YACYAXTROOT}"
echo "executables = ${executables}"
echo "### ------------------------------------------- ###"
### ---------------------------------------------------- ###

### --------------------- build CLEO ------------------- ###
buildcmd="${CLEO_PATH2BUILD}/scripts/bash/build_cleo.sh"
echo ${buildcmd}
eval ${buildcmd}
### ---------------------------------------------------- ###

### ---------------- compile executables --------------- ###
make_clean=true
compilecmd="${CLEO_PATH2BUILD}/scripts/bash/compile_cleo.sh ${executables} ${make_clean}"
echo ${compilecmd}
eval ${compilecmd}
### ---------------------------------------------------- ###
