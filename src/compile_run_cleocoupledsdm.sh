#!/bin/bash
#SBATCH --job-name=run_cleocoupledsdm
#SBATCH --partition=gpu
#SBATCH --gpus=4
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --mem=30G
#SBATCH --time=00:30:00
#SBATCH --mail-user=clara.bayley@mpimet.mpg.de
#SBATCH --mail-type=FAIL
#SBATCH --account=mh1126
#SBATCH --output=./gpu_cleocoupledsdm_out.%j.out
#SBATCH --error=./gpu_cleocoupledsdm_err.%j.out

### ------------- PLEASE NOTE: this script assumes you ------------- ###
### ------------- have already built CLEO in path2build ------------ ###
### -------------------  directory using cmake  -------------------- ###

### ----- You need to edit these lines to set your ----- ###
### ----- default compiler and python environment   ---- ###
### ----  and paths for CLEO and build directories  ---- ###
module load python3/2022.01-gcc-11.2.0

executable="cleocoupledsdm"
path2CLEO=${HOME}/CLEO/
path2build=$1              # get from command line argument
configfile=${HOME}/CLEO/src/config/config.txt
### ---------------------------------------------------- ###

### ----------------- compile executable --------------- ###
if [ "${path2build}" == "" ]
then
  path2build=${HOME}/CLEO/build/
fi

compilecmd="${path2CLEO}/scripts/bash/compile_cleo.sh ${path2build} ${executable}"
echo ${compilecmd}
${compilecmd}

### ------------------- run executable ----------------- ###
export OMP_PROC_BIND=spread
export OMP_PLACES=threads
runcmd="${path2build}/src/${executable} ${configfile}"
echo ${runcmd}
${runcmd}
### -------------------------------------------------- ###
