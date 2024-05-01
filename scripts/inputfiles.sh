#!/bin/bash
#SBATCH --job-name=inputfiles
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --mem=30G
#SBATCH --time=00:05:00
#SBATCH --mail-user=clara.bayley@mpimet.mpg.de
#SBATCH --mail-type=FAIL
#SBATCH --account=mh1126
#SBATCH --output=./inputfiles_out.%j.out
#SBATCH --error=./inputfiles_err.%j.out

### ----- You need to edit these lines to set your ----- ###
### ----- default compiler and python environment   ---- ###
### ----  and paths for CLEO and build directories  ---- ###
condaenv=/work/mh1126/m300950/cleoenv
module load python3/2022.01-gcc-11.2.0
source activate ${condaenv}

path2CLEO=${2:-${HOME}/CLEO}
path2scripts=${path2CLEO}/scripts
python=${condaenv}/bin/python

### ---------------------------------------------------- ###

configfile=$1
path2build=${3:-${path2CLEO}/build}

if [ "${configfile}" == "" ]
then
  echo "Please specify config file"
else

  echo "config file: ${configfile}"
  echo "path to build directory: ${path2build}"

  ### --------------- create gbx boundaries -------------- ###
  echo "${python} create_gbxboundariesbinary_script.py ${path2CLEO} ${path2build} ${configfile}"
  ${python} ${path2scripts}/create_gbxboundariesbinary_script.py ${path2CLEO} ${path2build} ${configfile}
  ### ---------------------------------------------------- ###

  ### -------- create superdrop initial conditions ------- ###
  echo "${python} create_initsuperdropsbinary_script.py ${path2CLEO} ${path2build} ${configfile}"
  ${python} ${path2scripts}/create_initsuperdropsbinary_script.py ${path2CLEO} ${path2build} ${configfile}
  ### ---------------------------------------------------- ###

  ### --------- create thermodynamics (optional) --------- ###
  echo "${python} create_thermobinaries_script.py ${path2CLEO} ${path2build} ${configfile}"
  ${python} ${path2scripts}/create_thermobinaries_script.py ${path2CLEO} ${path2build} ${configfile}
  ### ---------------------------------------------------- ###
fi
