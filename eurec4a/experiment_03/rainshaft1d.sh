#!/bin/bash
#SBATCH --job-name=sdm_eueurec4a_rain1d
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --mem=30G
#SBATCH --time=00:30:00
#SBATCH --mail-user=nils-ole.niebaum@mpimet.mpg.de
#SBATCH --mail-type=FAIL
#SBATCH --account=mh1126
#SBATCH --output=./rain1d_out.%j.out
#SBATCH --error=./rain1d_err.%j.out

### ----- You need to edit these lines to set your ----- ###
### ----- default compiler and python environment   ---- ###
### ----  and paths for CLEO and build directories  ---- ###
module load gcc/11.2.0-gcc-11.2.0
module load python3/2022.01-gcc-11.2.0
module load nvhpc/23.7-gcc-11.2.0
spack load cmake@3.23.1%gcc
source activate /work/mh1126/m300950/condaenvs/superdropsenv
logfile=${HOME}/CLEO/results/logging/rain1d_yaml.log

echo "============================================"
echo "START NEW RUN"
date
echo "git hash: $(git rev-parse HEAD)"
echo "--------------------------------------------"

path2CLEO=${HOME}/CLEO/
path2build=${HOME}/CLEO/build/
path2sdm_eurec4a=${HOME}/repositories/sdm-eurec4a/
rawdirectory=${path2CLEO}/data/output/raw/rain/
processeddirectory=${path2CLEO}/data/output/processed/rain/

configfile=${path2CLEO}eurec4a/experiment_03/src/config/rain1d_config.txt

yamldirectory=${HOME}/repositories/sdm-eurec4a/data/model/input/all_rain_clusters

python=/work/mh1126/m300950/condaenvs/superdropsenv/bin/python
pythonPySD=/work/mh1126/m301096/conda/envs/sdm_pysd_env312/bin/python
gxx="g++"
gcc="gcc"
cuda="nvc++"
### ---------------------------------------------------- ###

### ------------ choose Kokkos configuration ----------- ###
kokkosflags="-DKokkos_ARCH_NATIVE=ON -DKokkos_ARCH_AMPERE80=ON -DKokkos_ENABLE_SERIAL=ON" # serial kokkos
kokkoshost="-DKokkos_ENABLE_OPENMP=ON"                                          # flags for host parallelism (e.g. using OpenMP)
kokkosdevice="-DKokkos_ENABLE_CUDA=ON -DKokkos_ENABLE_CUDA_LAMBDA=ON"           # flags for device parallelism (e.g. using CUDA)
### ---------------------------------------------------- ###

# ### ------------------------ build --------------------- ###
### build CLEO using cmake (with openMP thread parallelism through Kokkos)
buildcmd="CXX=${gxx} CC=${gcc} CUDA=${cuda} cmake -S ${path2CLEO} -B ${path2build} ${kokkosflags} ${kokkoshost} ${kokkosdevice}"
echo ${buildcmd}
CXX=${gxx} CC=${gcc} CUDA=${cuda} cmake -S ${path2CLEO} -B ${path2build} ${kokkosflags} ${kokkoshost} ${kokkosdevice}

export OMP_PROC_BIND=spread
export OMP_PLACES=threads

### ensure these directories exist (it's a good idea for later use)
mkdir ${path2build}bin
mkdir ${path2build}share
### ---------------------------------------------------- ###

### ------------------- compile & run ------------------ ###
### generate input files and run 1-D rainshaft example

for yamlfile in ${yamldirectory}/*.yaml; do
    echo "Running rainshaft1d.py with ${yamlfile}"
    {
        ${python} rainshaft1d.py ${path2CLEO} ${path2build} ${configfile} ${yamlfile} ${rawdirectory} > ${logfile} 
    } || {
        echo "NO DATA CREATED"
        }
    # {
    #     ${pythonPySD} zarrfiles_to_netcdf.py ${path2CLEO} ${path2sdm_eurec4a} ${yamlfile} ${rawdirectory} ${processeddirectory}
    # } || {
    #     echo "NO DATA CREATED"
    # }
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
