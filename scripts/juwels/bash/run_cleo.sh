#!/bin/bash
#SBATCH --job-name=run_cleo
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=940M
#SBATCH --time=00:05:00
#SBATCH --mail-user=clara.bayley@mpimet.mpg.de
#SBATCH --mail-type=FAIL
#SBATCH --account=exaww
#SBATCH --output=./run_cleo_out.%j.out
#SBATCH --error=./run_cleo_err.%j.out

### Please note: script may assume required CLEO_[XXX]
### variables have already exported (!)

set -e
module purge

executable2run=$1
configfile=$2
stacksize_limit=$3 # kB
bashsrc=${CLEO_PATH2CLEO}/scripts/juwels/bash/src
### -------------------- check inputs ------------------ ###
source ${bashsrc}/check_inputs.sh
check_args_not_empty "${executable2run}" "${configfile}" "${CLEO_ENABLEYAC}" "${CLEO_COMPILERNAME}"
### ---------------------------------------------------- ###

### ----------- load compiler(s) and libraries --------- ###
source ${bashsrc}/juwels_packages.sh

if [ "${CLEO_COMPILERNAME}" == "intel" ]
then
  module load ${juwels_intel}
  module load ${juwels_intel_mpi}
elif [ "${CLEO_COMPILERNAME}" == "gcc" ]
then
  module load ${juwels_gcc}
  module load ${juwels_gcc_mpi}
  if [ "${CLEO_BUILDTYPE}" == "cuda" ]
  then
    echo "Bad inputs, CUDA build enabled but building CLEO with CUDA on JUWELS is not currently supported"
    exit 1
  fi
fi
### ---------------------------------------------------- ###

### ----------------- run executable --------------- ###
source ${bashsrc}/runtime_settings.sh ${stacksize_limit}
runcmd="${executable2run} ${configfile}"
echo ${runcmd}
eval ${runcmd}
### ---------------------------------------------------- ###
