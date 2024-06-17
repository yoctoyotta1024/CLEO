#!/bin/bash
#SBATCH --job-name=runexample
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --mem=30G
#SBATCH --time=00:10:00
#SBATCH --mail-user=clara.bayley@mpimet.mpg.de
#SBATCH --mail-type=FAIL
#SBATCH --account=mh1126
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
enableyac=$4      # required "true" or otherwise
executables="$5"
pythonscript=$6
script_args="$7"

cleoenv=/work/mh1126/m300950/cleoenv
python=${cleoenv}/bin/python3
yacyaxtroot=/work/mh1126/m300950/yac
spack load cmake@3.23.1%gcc
module load python3/2022.01-gcc-11.2.0
source activate ${cleoenv}
### ---------------------------------------------------- ###
### ---------------------------------------------------- ###
### ---------------------------------------------------- ###

### -------------------- print inputs ------------------ ###
echo "----- Running Example -----"
echo "buildtype:  ${buildtype}"
echo "path2CLEO: ${path2CLEO}"
echo "path2build: ${path2build}"
echo "enableyac: ${enableyac}"
echo "executables: ${executables}"
echo "pythonscript: ${pythonscript}"
echo "script_args: ${script_args}"
echo "---------------------------"
### ---------------------------------------------------- ###

### ---------------------- build CLEO ------------------ ###
${path2CLEO}/scripts/bash/build_cleo.sh ${buildtype} ${path2CLEO} ${path2build} ${enableyac} ${yacyaxtroot}
### ---------------------------------------------------- ###

### --------- compile executable(s) from scratch ---------- ###
cd ${path2build} && make clean
${path2CLEO}/scripts/bash/compile_cleo.sh ${cleoenv} ${buildtype} ${path2build} "${executables}"
### ---------------------------------------------------- ###

### --------- run model through Python script ---------- ###
export OMP_PROC_BIND=spread
export OMP_PLACES=threads

# TODO(all): add exports to paths required if YAC is enabled

${python} ${pythonscript} ${path2CLEO} ${path2build} ${script_args}
### ---------------------------------------------------- ###
