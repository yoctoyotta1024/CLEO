#!/bin/bash
#SBATCH --job-name=build_compile_cleo
#SBATCH --partition=compute
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

### ------------------ input parameters ---------------- ###
### ----- You need to edit these lines to specify ------ ###
### ----- (your environment and) directory paths ------- ###
### ---------- and executable(s) to compile ------------ ###
buildtype=$1                              # "serial", "threads", "openmp" or "cuda"
compilertype=${2:-intel}                  # "intel" or "gcc"
enableyac=${3:-false}                     # "true" or otherwise false
path2CLEO=${4:-${HOME}/CLEO}              # must be absolute path
path2build=${5:-${path2CLEO}/build}       # should be absolute path
executables=${6:-"cleocoupledsdm rain"}   # list of executables to compile
yacyaxtroot=/work/bm1183/m300950/yacyaxt  # used if enableyac == "true"
### ---------------------------------------------------- ###

if [[ "${buildtype}" == "" || "${compilertype}" == "" || "${enableyac}" == "" ||
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

export CLEO_BUILD=${buildtype}
export CLEO_COMPILER=${compilertype}
export CLEO_ENABLEYAC=${enableyac}
export CLEO_PATH2CLEO=${path2CLEO}
export CLEO_PATH2BUILD=${path2build}

if [ ${enableyac} == "true" ]
then
  export CLEO_YACYAXTROOT=${yacyaxtroot}
fi

### -------------------- print inputs ------------------- ###
echo "### --------------- Build Inputs -------------- ###"
echo "CLEO_BUILD = ${CLEO_BUILD}"
echo "CLEO_COMPILER = ${CLEO_COMPILER}"
echo "CLEO_ENABLEYAC = ${CLEO_ENABLEYAC}"
echo "CLEO_PATH2CLEO = ${CLEO_PATH2CLEO}"
echo "CLEO_PATH2BUILD = ${CLEO_PATH2BUILD}"
echo "CLEO_YACYAXTROOT = ${CLEO_YACYAXTROOT}"
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
