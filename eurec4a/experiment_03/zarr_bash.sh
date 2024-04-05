#!/bin/bash
#SBATCH --job-name=zarrfiles2netcdf_eurec4a_rain1d_subset
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --mem=30G
#SBATCH --time=01:30:00
#SBATCH --mail-user=nils-ole.niebaum@mpimet.mpg.de
#SBATCH --mail-type=FAIL
#SBATCH --account=mh1126
#SBATCH --output=./rain1d_out.%j.out
#SBATCH --error=./rain1d_err.%j.out

### ----- You need to edit these lines to set your ----- ###
### ----- default compiler and python environment   ---- ###
### ----  and paths for CLEO and build directories  ---- ###
pythonPySD=/work/mh1126/m301096/conda/envs/sdm_pysd_env312/bin/python
source activate ${pythonPySD}
logfile=${HOME}/CLEO/results/logging/rain1d_yaml.log

echo "============================================"
echo "START NEW RUN"
date
echo "git hash: $(git rev-parse HEAD)"
echo "--------------------------------------------"

path2CLEO=${HOME}/CLEO/
path2build=${HOME}/CLEO/build/
path2sdm_eurec4a=${HOME}/repositories/sdm-eurec4a/
path2experiment=${path2CLEO}/eurec4a/experiment_03/
rawdirectory=${path2CLEO}/data/output/raw/rain/
processeddirectory=${path2CLEO}/data/output/processed/rain/
# yamldirectory=${HOME}/repositories/sdm-eurec4a/data/model/input/all_rain_clusters
yamldirectory=${HOME}/repositories/sdm-eurec4a/data/model/input/subset


for yamlfile in ${yamldirectory}/*.yaml; do
    echo "Running zarrfiles_to_netcdf.py with ${yamlfile}"
    # {
    #     ${python} rainshaft1d.py ${path2CLEO} ${path2build} ${configfile} ${yamlfile} ${rawdirectory} > ${logfile} 
    # } || {
    #     echo "NO DATA CREATED"
    #     }
    {
        ${pythonPySD} ${path2experiment}/zarrfiles_to_netcdf.py ${path2CLEO} ${path2sdm_eurec4a} ${yamlfile} ${rawdirectory} ${processeddirectory}
    } || {
        echo "NO DATA CREATED"
    }
done
# echo "--------------------------------------------"
# echo "Plot results"
# source activate /work/mh1126/m301096/conda/envs/sdm_pysd_env312

# # ${pythonPySD} ${path2CLEO}eurec4a/experiment_02/plot_rainshaft1d.py ${yamlfile}

### ---------------------------------------------------- ###

echo "--------------------------------------------"
echo "END RUN"
date
echo "============================================"
