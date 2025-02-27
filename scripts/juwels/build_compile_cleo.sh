#!/bin/bash
#SBATCH --job-name=build_compile_cleo
#SBATCH --partition=devel
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=940M
#SBATCH --time=00:05:00
#SBATCH --mail-user=clara.bayley@mpimet.mpg.de
#SBATCH --mail-type=FAIL
#SBATCH --account=exaww
#SBATCH --output=./build_compile_cleo_out.%j.out
#SBATCH --error=./build_compile_cleo_err.%j.out

set -e
module purge

### ------------------ input parameters ---------------- ###
### ----- You need to edit these lines to specify ------ ###
### ----- your build configuration and executables ----- ###
### ---------------------------------------------------- ###
buildtype=$1                                # "serial", "threads", "openmp" or "cuda"
compilername=${2:-gcc}                      # "intel" or "gcc"
path2CLEO=${3:-${PROJECT}/bayley1/CLEO}     # must be absolute path
path2build=${4:-${path2CLEO}/build}         # should be absolute path
executables=${5:-"cleocoupledsdm"}          # list of executables to compile
enabledebug=${6:-false}                     # == "true" or otherwise false
enableyac=${7:-false}                       # == "true" or otherwise false
yacyaxtroot=${8:-NA}                        # yac and yaxt in yacyaxtroot/yac and yacyaxtroot/yaxt
enable_mptrac=${9:-true}                    # == "true" or otherwise false
make_clean=${10:-true}                      # == "true" or otherwise false
### ---------------------------------------------------- ###

### ------------------ check arguments ----------------- ###
if [ "${path2CLEO}" == "" ]
then
  echo "Please provide path to CLEO source directory"
  exit 1
fi
source ${path2CLEO}/scripts/juwels/bash/src/check_inputs.sh
check_args_not_empty "${buildtype}" "${compilername}" "${path2CLEO}" "${path2build}"
check_args_not_empty "${enabledebug}" "${enableyac}" "${enable_mptrac}"
### ---------------------------------------------------- ###

### ----------------- export inputs -------------------- ###
export CLEO_BUILDTYPE=${buildtype}
export CLEO_COMPILERNAME=${compilername}
export CLEO_PATH2CLEO=${path2CLEO}
export CLEO_PATH2BUILD=${path2build}
export CLEO_ENABLEDEBUG=${enabledebug}
export CLEO_ENABLEYAC=${enableyac}
export CLEO_ENABLE_MPTRAC=${enable_mptrac}

if [ ${CLEO_ENABLEYAC} == "true" ]
then
  export CLEO_YACYAXTROOT=${yacyaxtroot}
fi
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
echo "CLEO_ENABLEDEBUG = ${CLEO_ENABLEDEBUG}"
echo "CLEO_ENABLEYAC = ${CLEO_ENABLEYAC}"
echo "CLEO_YACYAXTROOT = ${CLEO_YACYAXTROOT}"
echo "CLEO_ENABLE_MPTRAC = ${CLEO_ENABLE_MPTRAC}"
echo "executables = ${executables}"
echo "### ------------------------------------------- ###"
### ---------------------------------------------------- ###

### --------------------- build CLEO ------------------- ###
buildcmd="${CLEO_PATH2CLEO}/scripts/juwels/bash/build_cleo.sh"
echo ${buildcmd}
eval ${buildcmd}
### ---------------------------------------------------- ###

### ---------------- compile executables --------------- ###
compilecmd="${CLEO_PATH2CLEO}/scripts/juwels/bash/compile_cleo.sh \"${executables}\" ${make_clean}"
echo ${compilecmd}
eval ${compilecmd}
### ---------------------------------------------------- ###
