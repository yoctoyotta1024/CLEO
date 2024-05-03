#!/bin/bash
#SBATCH --job-name=eurec4a1d
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --gpus=4
#SBATCH --ntasks-per-node=128
#SBATCH --mem=30G
#SBATCH --time=00:10:00
#SBATCH --mail-user=nils-ole.niebaumy@mpimet.mpg.de
#SBATCH --mail-type=FAIL
#SBATCH --account=mh1126
#SBATCH --output=./logfiles/eurec4a1d_out.%j.out
#SBATCH --error=./logfiles/eurec4a1d_err.%j.out

### ---------------------------------------------------- ###
### ------------------ Input Parameters ---------------- ###
### ------ You MUST edit these lines to set your ------- ###
### ----- environment, build type, directories, the ---- ###
### --------- executable(s) to compile and your -------- ###
### --------------  python script to run. -------------- ###
### ---------------------------------------------------- ###
buildtype="cuda"
path2CLEO=${HOME}/CLEO/
path2build=${HOME}/CLEO/build_eurec4a1D/
executables="eurec4a1D"
pythonscript=${path2CLEO}examples/eurec4a1d/eurec4a1d.py

configfile=${path2CLEO}examples/eurec4a1d/src/config/eurec4a1d_config.yaml
rawdirectory=${path2CLEO}data/output/raw/no_aerosols_no_collision/

path2sdmeurec4a=${HOME}/repositories/sdm-eurec4a/
cloud_observation_configfile=${path2sdmeurec4a}data/model/input/new/clusters_18.yaml


script_args="${HOME} ${configfile} ${cloud_observation_configfile} ${rawdirectory}"


cleoenv=/work/mh1126/m300950/cleoenv
python=${cleoenv}/bin/python3
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
echo "executables: ${executables}"
echo "pythonscript: ${pythonscript}"
echo "script_args: ${script_args}"
echo "---------------------------"
### ---------------------------------------------------- ###

### ---------------------- build CLEO ------------------ ###
${path2CLEO}/scripts/bash/build_cleo.sh ${buildtype} ${path2CLEO} ${path2build}
### ---------------------------------------------------- ###

### --------- compile executable(s) from scratch ---------- ###
cd ${path2build} && make clean

${path2CLEO}/scripts/bash/compile_cleo.sh ${cleoenv} ${buildtype} ${path2build} "${executables}"
### ---------------------------------------------------- ###

### --------- run model through Python script ---------- ###
export OMP_PROC_BIND=spread
export OMP_PLACES=threads
${python}  ${pythonscript} ${path2CLEO} ${path2build} ${script_args}
### ---------------------------------------------------- ###

echo "--------------------------------------------"
echo "END RUN"
date
echo "============================================"
