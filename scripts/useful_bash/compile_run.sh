#!/bin/bash
#SBATCH --job-name=runCLEO
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --mem=30G
#SBATCH --time=00:30:00
#SBATCH --mail-user=clara.bayley@mpimet.mpg.de
#SBATCH --mail-type=FAIL
#SBATCH --account=mh1126
#SBATCH --output=./build/bin/runCLEO_out.%j.out
#SBATCH --error=./build/bin/runCLEO_err.%j.out

### ------------- PLEASE NOTE: this script assumes you ------------- ###
### ------------- have already built CLEO in path2build ------------ ### 
### -------------------  directory using cmake  -------------------- ###

### ----- You need to edit these lines to set your ----- ###
### ----- default compiler and python environment   ---- ###
### ----  and paths for CLEO and build directories  ---- ###
module load gcc/11.2.0-gcc-11.2.0
module load python3/2022.01-gcc-11.2.0
source activate /work/mh1126/m300950/condaenvs/cleoenv 
path2CLEO=${HOME}/CLEO/
path2build=${HOME}/CLEO/build/
configfile=${HOME}/CLEO/src/config/config.txt
python=python

# path2CLEO=${HOME}/Documents/b2_springsummer2023/CLEO/
# path2build=${HOME}/Documents/b2_springsummer2023/CLEO/build/
# python=${HOME}/opt/anaconda3/envs/superdropsV2/bin/python
### ---------------------------------------------------- ###

### ------------------- compile_run.sh ----------------- ###
### ensure these directories exist (it's a good idea for later use)
mkdir ${path2build}bin
mkdir ${path2build}share

### compile CLEO in ./build directory
cd ${path2build} && pwd 
make -j 128

# ### generate input files
# ${python} ${path2CLEO}create_gbxboundariesbinary_script.py ${path2CLEO} ${path2build} ${configfile}
# ${python} ${path2CLEO}create_thermobinaries_script.py ${path2CLEO} ${path2build} ${configfile}
# ${python} ${path2CLEO}create_initsuperdropsbinary_script.py ${path2CLEO} ${path2build} ${configfile}

### run CLEO
export OMP_PROC_BIND=spread
export OMP_PLACES=threads
runcmd="${path2build}/src/runCLEO ${configfile}"
echo ${runcmd}
${runcmd}
### ---------------------------------------------------- ###