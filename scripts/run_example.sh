#!/bin/bash
#SBATCH --job-name=runexample
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --mem=30G
#SBATCH --time=00:30:00
#SBATCH --mail-user=clara.bayley@mpimet.mpg.de
#SBATCH --mail-type=FAIL
#SBATCH --account=mh1126
#SBATCH --output=./runexample_out.%j.out
#SBATCH --error=./runexample_err.%j.out

example=$1
sbatch=$2
examplesdir=${3:-./examples}

if [ "${example}" == "" ]
then
  echo "Please specify and example to run"
  exit 0

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

  elif  [ "${example}" == "constthermo2d" ]
    then
    ${sbatch} ${examplesdir}/constthermo2d/constthermo2d.sh

  elif  [ "${example}" == "divfree2d" ]
  then
    ${sbatch} ${examplesdir}/divfreemotion/divfree2d.sh

  elif  [ "${example}" == "eurec4a1d" ]
  then
    ${sbatch} ${examplesdir}/examples/eurec4a1d/eurec4a1d.sh

  elif  [ "${example}" == "rainshaft1d" ]
  then
    ${sbatch} ${examplesdir}/examples/rainshaft1d/rainshaft1d.sh

  elif  [ "${example}" == "speedtest" ]
  then
    ${sbatch} ${examplesdir}/ examples/speedtest/speedtest.sh

  else
    echo "'${example}' is not an example"
  fi
fi
