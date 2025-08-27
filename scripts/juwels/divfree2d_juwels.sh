#!/bin/bash
#SBATCH --job-name=divfree2d
#SBATCH --partition=devel
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=10G
#SBATCH --time=00:10:00
#SBATCH --mail-user=clara.bayley@mpimet.mpg.de
#SBATCH --mail-type=FAIL
#SBATCH --account=exaww
#SBATCH --output=./divfree2d_out.%j.out
#SBATCH --error=./divfree2d_err.%j.out

### ------------------ Input Parameters ---------------- ###
### ------ You MUST edit these lines to set your ------- ###
### ---- build type, directories, the executable(s) ---- ###
### -------- to compile, and your python script -------- ###
### ---------------------------------------------------- ###
compilername=intel
buildtype="openmp"
path2CLEO=${PROJECT}/bayley1/CLEO/
path2build=${PROJECT}/bayley1/CLEO/build_divfree2d/
build_flags="-DCLEO_COUPLED_DYNAMICS=fromfile -DCLEO_DOMAIN=cartesian -DCLEO_NO_PYBINDINGS=true"
executables="divfree2d"

python=/home/m/m300950/CLEO/.venv/bin/python3
pythonscript=${path2CLEO}/examples/divfreemotion/divfree2d.py
configfile=${path2CLEO}/examples/divfreemotion/src/config/divfree2d_config.yaml
script_args="${configfile}"

yacyaxtroot=NA
enabledebug=false
make_clean=false
stacksize_limit=204800 # ulimit -s [stacksize_limit] (kB)
### ---------------------------------------------------- ###
### ---------------------------------------------------- ###
### ---------------------------------------------------- ###

### -------------------- print inputs ------------------ ###
echo "----- Running Example -----"
echo "buildtype = ${buildtype}"
echo "compilername = ${compilername}"
echo "path2CLEO = ${path2CLEO}"
echo "path2build = ${path2build}"
echo "yacyaxtroot = ${yacyaxtroot}"
echo "executables = ${executables}"
echo "pythonscript = ${pythonscript}"
echo "script_args = ${script_args}"
echo "---------------------------"
### ---------------------------------------------------- ###

### --------------- build and compile CLEO ------------- ###
cmd="${path2CLEO}/scripts/juwels/build_compile_cleo.sh \
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
### ---------------------------------------------------- ###

### ----------- load compiler(s) and libraries --------- ###
source ${bashsrc}/juwels_packages.sh

if [ "${compilername}" == "intel" ]
then
  module load ${juwels_intel}
  module load ${juwels_intel_mpi}
elif [ "${compilername}" == "gcc" ]
then
  module load ${juwels_gcc}
  module load ${juwels_gcc_mpi}
  if [ "${CLEO_BUILDTYPE}" == "cuda" ]
  then
    echo "Bad inputs, CUDA build enabled but building CLEO with CUDA on JUWELS is not currently supported"
    exit 1
  fi
fi
### ---------------------------------------------------- ###

### --------- run model through Python script ---------- ###
export CLEO_PATH2CLEO=${path2CLEO}
export CLEO_BUILDTYPE=${buildtype}
export CLEO_COMPILERNAME=${compilername}
export CLEO_YACYAXTROOT=${yacyaxtroot}
source ${path2CLEO}/scripts/juwels/bash/src/runtime_settings.sh ${stacksize_limit}

# TODO(ALL): split python scripts away from running executable
${python} ${pythonscript} ${path2CLEO} ${path2build} ${script_args}
### ---------------------------------------------------- ###
