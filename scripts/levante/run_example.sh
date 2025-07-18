#!/bin/bash

example=$1
sbatch=$2
path2CLEO=${3:-${HOME}/CLEO}
examplesdir=${path2CLEO}/examples

if [ "${example}" == "" ]
then
  echo "Please specify and example to run"
  exit 1

else
  if [ "${example}" == "as2017" ]
  then
    ${sbatch} ${examplesdir}/adiabaticparcel/as2017.sh

  elif  [ "${example}" == "cuspbifurc" ]
  then
    ${sbatch} ${examplesdir}/adiabaticparcel/cuspbifurc.sh

  elif  [ "${example}" == "shima2009" ]
    then
    ${sbatch} ${examplesdir}/boxmodelcollisions/shima2009.sh

  elif  [ "${example}" == "breakup" ]
    then
    ${sbatch} ${examplesdir}/boxmodelcollisions/breakup.sh

  elif  [ "${example}" == "bubble3d" ]
    then
      ${sbatch} ${examplesdir}/bubble3d/bubble3d.sh # WIP

  elif  [ "${example}" == "constthermo2d" ]
    then
    ${sbatch} ${examplesdir}/constthermo2d/constthermo2d.sh

  elif  [ "${example}" == "divfree2d" ]
  then
    ${sbatch} ${examplesdir}/divfreemotion/divfree2d.sh

  elif  [ "${example}" == "eurec4a1d" ]
  then
    ${sbatch} ${examplesdir}/eurec4a1d/eurec4a1d.sh

  elif  [ "${example}" == "python_bindings" ]
    then
      ${sbatch} ${examplesdir}/python_bindings/python_bindings.sh

  elif  [ "${example}" == "rainshaft1d" ]
  then
    ${sbatch} ${examplesdir}/rainshaft1d/rainshaft1d.sh

  elif  [ "${example}" == "speedtest" ]
  then
    ${sbatch} ${examplesdir}/speedtest/speedtest.sh

  else
    echo "'${example}' is not an example"
    exit 1
  fi
fi
