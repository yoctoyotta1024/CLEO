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
#SBATCH --output=./run_cleocoupledsdm_out.%j.out
#SBATCH --error=./run_cleocoupledsdm_err.%j.out

### ------------- PLEASE NOTE: this script assumes you have ------------- ###
### ------- already built CLEO in "path2build" and compiled the --------- ###
### ------------- desired executable in "path2build/src" ---------------- ###

### ------------------ input parameters ---------------- ###
### ----- You need to edit these lines to specify ------ ###
### ------ (your environment and) directory paths ------ ###
executable="cleocoupledsdm"
path2CLEO=${HOME}/CLEO/
path2build=$1 # get from command line argument
configfile=${HOME}/CLEO/src/config/config.txt

if [ "${path2build}" == "" ]
then
  path2build=${HOME}/CLEO/build/
fi
### ---------------------------------------------------- ###

### ----------------- compile executable --------------- ###
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
