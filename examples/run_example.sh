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
buildtype=$1
path2CLEO=$2
path2build=$3
enableyac=$4
executables="$5"
pythonscript=$6
script_args="$7"

cleoenv=/work/bm1183/m300950/bin/envs/cleoenv
python=${cleoenv}/bin/python3
enabledebug=false
make_clean=false
yacyaxtroot=/work/bm1183/m300950/yacyaxt
stacksize_limit=204800 # ulimit -s [stacksize_limit] (kB)

if [ "${buildtype}" == "cuda" ]
then
  compilername=gcc
else
  compilername=intel
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
echo "enableyac = ${enableyac}"
echo "executables = ${executables}"
echo "pythonscript = ${pythonscript}"
echo "script_args = ${script_args}"
echo "---------------------------"
### ---------------------------------------------------- ###

### --------------- build and compile CLEO ------------- ###
cmd="${path2CLEO}/scripts/build_compile_cleo.sh \
  ${buildtype}
  ${compilername}
  ${path2CLEO}
  ${path2build}
  "\"${executables}\""
  ${enabledebug}
  ${enableyac}
  ${yacyaxtroot}
  ${make_clean}"
echo ${cmd}
eval ${cmd}
### ---------------------------------------------------- ###

### --------- run model through Python script ---------- ###
export CLEO_PATH2CLEO=${path2CLEO}
export CLEO_BUILDTYPE=${buildtype}
export CLEO_ENABLEYAC=${enableyac}
source ${path2CLEO}/scripts/bash/src/runtime_settings.sh ${stacksize_limit}

# TODO(all): split python scripts away from running executable
${python} ${pythonscript} ${path2CLEO} ${path2build} ${script_args}
### ---------------------------------------------------- ###
