#!/bin/bash
#SBATCH --job-name=speedtest
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --gpus=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=256
#SBATCH --mem=10G
#SBATCH --time=00:30:00
#SBATCH --mail-user=clara.bayley@mpimet.mpg.de
#SBATCH --mail-type=FAIL
#SBATCH --account=bm1183
#SBATCH --output=./speedtest_out.%j.out
#SBATCH --error=./speedtest_err.%j.out

### ------------------ Input Parameters ---------------- ###
### ------ You MUST edit these lines to set your ------- ###
### ---- build type, directories, the executable(s) ---- ###
### -------- to compile, and your python script -------- ###
### ---------------------------------------------------- ###
path2CLEO=${HOME}/CLEO/
path2build=${HOME}/CLEO/build_spdtest/
path2kokkostools=/work/bm1183/m300950/kokkos_tools_lib/lib64/
enableyac=false
executables="spdtest"

pythonscript=${path2CLEO}/examples/speedtest/speedtest.py
configfile=${path2CLEO}/examples/speedtest/src/config/speedtest_config.yaml
### ---------------------------------------------------- ###
### ---------------------------------------------------- ###
### ---------------------------------------------------- ###

# ensure these directories exist (it's a good idea for later use)
mkdir ${path2build}
mkdir ${path2build}/bin

### ---- run test for different types of parallelism ---- ###
buildtypes=("cuda" "openmp" "serial")
for buildtype in "${buildtypes[@]}"
do
  script_args="${configfile} ${path2build}/bin/ ${path2kokkostools} ${buildtype}"
  path2build_test=${path2build}${buildtype}"/"

  echo "build type: ${buildtype}"
  echo "path2build: ${path2build_test}"

  mkdir ${path2build_test}/bin
  mkdir ${path2build_test}/share

  ### ---------- build, compile and run example ---------- ###
  ${path2CLEO}/examples/run_example.sh \
    ${buildtype} ${path2CLEO} ${path2build_test} ${enableyac} \
    "${executables}" ${pythonscript} "${script_args}"
  ### ---------------------------------------------------- ###
done
### ---------------------------------------------------- ###
