#!/bin/bash

### ----- You need to edit these lines to set your ----- ###
### ----- default compiler and python environment   ---- ###
### ----  and paths for CLEO and build directories  ---- ###
configfile=$1
path2CLEO=${2:-${HOME}/CLEO}
path2build=${3:-${path2CLEO}/build}

path2scripts=${path2CLEO}/scripts
python=/home/m/m300950/CLEO/.venv/bin/python3
### ---------------------------------------------------- ###

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
