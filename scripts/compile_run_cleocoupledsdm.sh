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

### ------------- PLEASE NOTE: this script assumes you ------------- ###
### ----------- have already built CLEO in "path2build" ------------ ###
### -------------------  directory using cmake  -------------------- ###

### ------------------ input parameters ---------------- ###
### ----- You need to edit these lines to specify ------ ###
### ----- (your environment and) directory paths ------- ###
### ------------ and executable to compile ------------- ###

cleoenv=/work/mh1126/m300950/cleoenv
buildtype=$1
path2CLEO=${2:-${HOME}/CLEO}
path2build=${3:-${path2CLEO}/build} # get from command line argument
executable="cleocoupledsdm"
configfile=${path2CLEO}/roughpaper/src/config/config.yaml
run_executable=${path2build}/roughpaper/src/${executable}

### ---------------------------------------------------- ###

### ----------------- compile executable --------------- ###
rm ${run_executable}
compilecmd="${path2CLEO}/scripts/bash/compile_cleo.sh ${cleoenv} ${buildtype} ${path2build} ${executable}"
echo ${compilecmd}
${compilecmd}
### ---------------------------------------------------- ###

### ------------------- run executable ----------------- ###
cd ${path2build} && pwd
runcmd="${path2CLEO}/scripts/bash/run_cleo.sh ${run_executable} ${configfile}"
echo ${runcmd}
${runcmd}
### -------------------------------------------------- ###
