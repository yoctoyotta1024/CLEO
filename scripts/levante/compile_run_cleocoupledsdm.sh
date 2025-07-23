#!/bin/bash
#SBATCH --job-name=compile_run_cleo
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --gpus=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=940M
#SBATCH --time=00:30:00
#SBATCH --mail-user=clara.bayley@mpimet.mpg.de
#SBATCH --mail-type=FAIL
#SBATCH --account=bm1183
#SBATCH --output=./compile_run_cleo_out.%j.out
#SBATCH --error=./compile_run_cleo_err.%j.out

set -e
source /etc/profile
module purge
spack unload --all

### ------------------ input parameters ---------------- ###
### ----- You need to edit these lines to specify ------ ###
### ----- your build configuration and executables ----- ###
### ---------------------------------------------------- ###
buildtype=$1                                                     # "serial", "threads", "openmp" or "cuda"
compilername=${2:-intel}                                         # "intel" or "gcc"
path2CLEO=${3:-${HOME}/CLEO}                                     # must be absolute path
path2build=${4:-${path2CLEO}/build}                              # should be absolute path
yacyaxtroot=${5:-/work/bm1183/m300950/yacyaxt/${compilername}}   # yac and yaxt in yacyaxtroot/yac and yacyaxtroot/yaxt
executables=${6:-"cleocoupledsdm"}                               # executable(s) to compile or "NONE"
executable2run=${7:-${path2build}/roughpaper/src/${executables}} # path to executable to run
configfile=${8:-${path2CLEO}/roughpaper/src/config/config.yaml}  # configuration to run
stacksize_limit=${9:-204800}                                     # ulimit -s [stacksize_limit] (kB)
### ---------------------------------------------------- ###

### -------------------- check inputs ------------------ ###
if [[ "${buildtype}" == "" || "${compilername}" == "" ||
      "${path2CLEO}" == "" || "${path2build}" == "" ]]
then
  echo "Bad inputs, please check all the required inputs have been specified"
  exit 1
fi

if [[ "${path2CLEO}" == "${path2build}" ]]
then
  echo "Bad inputs, build directory cannot match the path to CLEO source"
  exit 1
fi

if [ "${buildtype}" != "serial" ] &&
   [ "${buildtype}" != "openmp" ] &&
   [ "${buildtype}" != "threads" ] &&
   [ "${buildtype}" != "cuda" ];
then
  echo "Bad inputs, build type must be 'serial', 'openmp', 'threads' or 'cuda'"
  exit 1
fi
### ---------------------------------------------------- ###

### ----------------- export inputs -------------------- ###
export CLEO_BUILDTYPE=${buildtype}
export CLEO_COMPILERNAME=${compilername}
export CLEO_PATH2CLEO=${path2CLEO}
export CLEO_PATH2BUILD=${path2build}
export CLEO_YACYAXTROOT=${yacyaxtroot}
### ---------------------------------------------------- ###

### --------------- print compiling inputs ------------- ###
echo "### --------------- User Inputs -------------- ###"
echo "CLEO_BUILDTYPE = ${CLEO_BUILDTYPE}"
echo "CLEO_COMPILERNAME = ${CLEO_COMPILERNAME}"
echo "CLEO_PATH2BUILD = ${CLEO_PATH2BUILD}"
echo "CLEO_YACYAXTROOT = ${CLEO_YACYAXTROOT}"
echo "### ------------------------------------------- ###"
### ---------------------------------------------------- ###

### ---------------- compile executables --------------- ###
make_clean=false
rm -f ${executable2run}
compilecmd="${CLEO_PATH2CLEO}/scripts/levante/bash/compile_cleo.sh ${executables} ${make_clean}"
echo ${compilecmd}
eval ${compilecmd}
### ---------------------------------------------------- ###

### -------------- print running inputs ---------------- ###
echo "### --------------- User Inputs -------------- ###"
echo "CLEO_COMPILERNAME = ${CLEO_COMPILERNAME}"
echo "CLEO_YACYAXTROOT = ${CLEO_YACYAXTROOT}"
echo "executable = ${executable2run}"
echo "config file for executable = ${configfile}"
echo "### ------------------------------------------- ###"
### ---------------------------------------------------- ###

### ------------------- run executable ----------------- ###
cd ${CLEO_PATH2BUILD} && pwd
runcmd="${CLEO_PATH2CLEO}/scripts/levante/bash/run_cleo.sh ${executable2run} ${configfile} ${stacksize_limit}"
echo ${runcmd}
eval ${runcmd}
### -------------------------------------------------- ###
