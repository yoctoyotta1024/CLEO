#!/bin/bash
#SBATCH --job-name=runexample
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=940MB
#SBATCH --time=00:10:00
#SBATCH --mail-user=clara.bayley@mpimet.mpg.de
#SBATCH --mail-type=FAIL
#SBATCH --account=bm1183
#SBATCH --output=./runexample_out.%j.out
#SBATCH --error=./runexample_err.%j.out

### ------ Generic script to build CLEO, compile  ------ ###
### ----- some of its executables and run a python ----- ###
### ------  script e.g. for example(s) on Levante. ----- ###

### ---------------------------------------------------- ###
### ------------------ Input Parameters ---------------- ###
### ------ You MUST edit these lines to set your ------- ###
### ----- environment, build type, directories, the ---- ###
### --------- executable(s) to compile and your -------- ###
### --------------  python script to run. -------------- ###
### ---------------------------------------------------- ###
do_build=$1  # == "true" or otherwise false
buildtype=$2
compilername=$3
path2CLEO=$4
path2build=$5
build_flags=$6
executables="$7"
pythonscript=$8
script_args="$9"

python=/home/m/m300950/CLEO/.venv/bin/python3
enabledebug=false
make_clean=false
stacksize_limit=204800 # ulimit -s [stacksize_limit] (kB)

if [[ "${buildtype}" == "cuda" && "${compilername}" != "gcc" ]];
then
  echo "CUDA build on Levante currently only compatible with gcc compiler"
  echo "-> please use compilername=gcc"
  exit 1
fi

if [[ "${compilername}" == "gcc" ]]
then
  yacyaxtroot=/work/bm1183/m300950/yacyaxt/gcc
elif [[ "${compilername}" == "intel" ]]
then
  yacyaxtroot=/work/bm1183/m300950/yacyaxt/intel
fi
### ---------------------------------------------------- ###
### ---------------------------------------------------- ###
### ---------------------------------------------------- ###

### -------------------- print inputs ------------------ ###
echo "----- Running Example -----"
echo "buildtype = ${buildtype}"
echo "compilername = ${compilername}"
echo "path2CLEO = ${path2CLEO}"
echo "path2build = ${path2build}"
echo "build_flags = ${build_flags}"
echo "yacyaxtroot = ${yacyaxtroot}"
echo "executables = ${executables}"
echo "pythonscript = ${pythonscript}"
echo "script_args = ${script_args}"
echo "---------------------------"
### ---------------------------------------------------- ###

### --------------- build and compile CLEO ------------- ###
if [ "${do_build}" == "true" ]
then
  cmd="${path2CLEO}/scripts/levante/build_compile_cleo.sh \
    ${buildtype}
    ${compilername}
    ${path2CLEO}
    ${path2build}
    ${yacyaxtroot}
    "\"${build_flags}\""
    "\"${executables}\""
    ${enabledebug}
    ${make_clean}"
  echo ${cmd}
  eval ${cmd}
fi
### ---------------------------------------------------- ###

### --------- run model through Python script ---------- ###
export CLEO_PATH2CLEO=${path2CLEO}
export CLEO_BUILDTYPE=${buildtype}
export CLEO_COMPILERNAME=${compilername}
export CLEO_YACYAXTROOT=${yacyaxtroot}
source ${path2CLEO}/scripts/levante/bash/src/runtime_settings.sh ${stacksize_limit}

# TODO(ALL): split python scripts away from running executable
${python} ${pythonscript} ${path2CLEO} ${path2build} ${script_args}
### ---------------------------------------------------- ###
