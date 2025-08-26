#!/bin/bash

example=$1
sbatch=$2
path2CLEO=${3:-${HOME}/CLEO}
path2examplesbash=${path2CLEO}/scripts/levante/examples

if [ "${example}" == "" ]
then
  echo "Please specify and example to run"
  exit 1

else
  if [ "${example}" == "as2017" ]
  then
    ${sbatch} ${path2examplesbash}/as2017.sh

  elif  [ "${example}" == "cuspbifurc" ]
  then
    ${sbatch} ${path2examplesbash}/cuspbifurc.sh

  elif  [ "${example}" == "shima2009" ]
    then
    ${sbatch} ${path2examplesbash}/shima2009.sh

  elif  [ "${example}" == "breakup" ]
    then
    ${sbatch} ${path2examplesbash}/breakup.sh

  elif  [ "${example}" == "bubble3d" ]
    then
      ${sbatch} ${path2examplesbash}/bubble3d.sh

  elif  [ "${example}" == "constthermo2d" ]
    then
    ${sbatch} ${path2examplesbash}/constthermo2d.sh

  elif  [ "${example}" == "divfree2d" ]
  then
    ${sbatch} ${path2examplesbash}/divfree2d.sh

  elif  [ "${example}" == "eurec4a1d" ]
  then
    ${sbatch} ${path2examplesbash}/eurec4a1d.sh

  elif  [ "${example}" == "fromfile" ]
  then
    ${sbatch} ${path2examplesbash}/fromfile.sh

  elif  [ "${example}" == "fromfile_irreg" ]
  then
    ${sbatch} ${path2examplesbash}/fromfile_irreg.sh

  elif  [ "${example}" == "python_bindings" ]
    then
      ${sbatch} ${path2examplesbash}/python_bindings.sh

  elif  [ "${example}" == "rainshaft1d" ]
  then
    ${sbatch} ${path2examplesbash}/rainshaft1d.sh

  elif  [ "${example}" == "kokkostools" ]
  then
    ${sbatch} ${path2examplesbash}/kokkostools.sh

  else
    echo "'${example}' is not an example"
    exit 1
  fi
fi
