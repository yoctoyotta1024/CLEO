#!/bin/bash
#SBATCH --job-name=build_cleocoupledsdm
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --gpus=4
#SBATCH --ntasks-per-node=128
#SBATCH --mem=30G
#SBATCH --time=00:05:00
#SBATCH --mail-user=clara.bayley@mpimet.mpg.de
#SBATCH --mail-type=FAIL
#SBATCH --account=mh1126
#SBATCH --output=./build/bin/build_cleocoupledsdm_out.%j.out
#SBATCH --error=./build/bin/build_cleocoupledsdm_err.%j.out

### ------------------ input parameters ---------------- ###
### ----- You need to edit these lines to specify ------ ###
### ----- (your environment and) directory paths ------- ###
### ---------- and executable(s) to compile ------------ ###
spack load cmake@3.23.1%gcc
cleoenv=/work/mh1126/m300950/cleoenv

buildtype=$1
enableyac=${2:-false}                     # "true" or otherwise
path2CLEO=${3:-${HOME}/CLEO}
path2build=${4:-${path2CLEO}/build}
yacyaxtroot=/work/mh1126/m300950/yac      # used if enableyac == "true"
executables="cleocoupledsdm"
### ---------------------------------------------------- ###

if [[ "${buildtype}" != "" && "${path2CLEO}" != "" && "${path2build}" != "" &&
    "${executables}" != "" && "${path2CLEO}" != "${path2build}" ]]
then

  if ! [ ${enableyac} == "true" ]
  then
    yacyaxtroot="" # don't provide path to YAC if build shouldn't require it
  fi

  ### --------------------- build CLEO ------------------- ###
  buildcmd="${path2CLEO}/scripts/bash/build_cleo.sh ${buildtype} ${path2CLEO} ${path2build} ${enableyac} ${yacyaxtroot}"
  echo ${buildcmd}
  ${buildcmd}
  ### ---------------------------------------------------- ###

  ### ---------------- compile executables --------------- ###
  cd ${path2build} && make clean
  compilecmd="${path2CLEO}/scripts/bash/compile_cleo.sh ${cleoenv} ${buildtype} ${path2build} ${executables}"
  echo ${compilecmd}
  ${compilecmd}
  ### ---------------------------------------------------- ###
else
  echo "Bad inputs, please check your buildtype, path2CLEO, path2build and executables names"
fi
