#!/bin/bash
#SBATCH --job-name=run_cleo
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --gpus=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=940M
#SBATCH --time=00:30:00
#SBATCH --account=bm1183
#SBATCH --output=./run_cleo_out.%j.out
#SBATCH --error=./run_cleo_err.%j.out

### Please note: script may assume required CLEO_[XXX]
### variables have already exported (!)

set -e
source /etc/profile
module purge
spack unload --all

executable2run=$1
configfile=$2
stacksize_limit=$3 # kB

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
bashsrc=${SCRIPT_DIR}/src

### -------------------- check inputs ------------------ ###
source ${bashsrc}/check_inputs.sh
check_args_not_empty "${executable2run}" "${configfile}"
check_args_not_empty "${CLEO_COMPILERNAME}" "${CLEO_YACYAXTROOT}"
### ---------------------------------------------------- ###

### ----------------- run executable --------------- ###
source ${bashsrc}/runtime_settings.sh ${stacksize_limit}
runcmd="${executable2run} ${configfile}"
echo ${runcmd}
eval ${runcmd}
### ---------------------------------------------------- ###
