#!/bin/bash
#SBATCH --job-name=eurec4a1d
#SBATCH --partition=compute
#SBATCH --nodes=1
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

echo "--------------------------------------------"
echo "START RUN"
date
echo "git hash: $(git rev-parse HEAD)"
echo "git branch: $(git symbolic-ref --short HEAD)"
echo "============================================"

buildtype="openmp"
path2CLEO=${HOME}/CLEO/
path2builds=${path2CLEO}builds/
path2data=${path2CLEO}data/output/raw/
path2eurec4a1d=${path2CLEO}examples/eurec4a1d/
# Use the stationary or evolving version of the model

### ---------- Setup for the EUREC4A1D model ---------- ###

# --- stationary version, with super droplet creation at domain top by boundarz conditions
path2build=${path2builds}build_eurec4a1D_stationary/
pythonscript=${path2eurec4a1d}scripts/eurec4a1d_stationary.py
configfile=${path2eurec4a1d}src/config/eurec4a1d_config_stationary.yaml
rawdirectory=${path2data}stationary/

# # --- evolving version, without super droplet creation at domain top
# path2build=${path2builds}build_eurec4a1D_evolving/
# pythonscript=${path2eurec4a1d}scripts/eurec4a1d_evolving.py
# configfile=${path2eurec4a1d}src/config/eurec4a1d_config_evolving.yaml
# rawdirectory=${path2data}evolving/


executables="eurec4a1D"


path2sdmeurec4a=${HOME}/repositories/sdm-eurec4a/
cloud_observation_configfile=${path2sdmeurec4a}data/model/input/new/clusters_18.yaml

# create the script arguments
script_args="${HOME} ${configfile} ${cloud_observation_configfile} ${rawdirectory}"
### ---------------------------------------------------- ###



### ------------------ Load Modules -------------------- ###
cleoenv=/work/mh1126/m300950/cleoenv
python=${cleoenv}/bin/python3
spack load cmake@3.23.1%gcc
module load python3/2022.01-gcc-11.2.0
source activate ${cleoenv}
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

# ### ---------------------- build CLEO ------------------ ###
# ${path2CLEO}/scripts/bash/build_cleo.sh ${buildtype} ${path2CLEO} ${path2build}
# ### ---------------------------------------------------- ###

# ### --------- compile executable(s) from scratch ---------- ###
# cd ${path2build} && make clean

# ${path2CLEO}/scripts/bash/compile_cleo.sh ${cleoenv} ${buildtype} ${path2build} "${executables}"
# ### ---------------------------------------------------- ###

### --------- run model through Python script ---------- ###
export OMP_PROC_BIND=spread
export OMP_PLACES=threads
${python}  ${pythonscript} ${path2CLEO} ${path2build} ${script_args}
### ---------------------------------------------------- ###

echo "--------------------------------------------"
echo "END RUN"
date
echo "============================================"
