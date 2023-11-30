#!/bin/bash

configfile=${HOME}/CLEO/src/config/config.txt

if [[ $1 == "1" ]];
then
  cd ${HOME}/CLEO/buildserial/ && pwd 
  runcmd="${HOME}/CLEO/buildserial/src/runCLEO ${configfile}"
  echo ${runcmd}
  ${runcmd}
  echo "serial run"
fi

if [[ $1 == "2" ]];
then

  export OMP_PROC_BIND=spread
  export OMP_PLACES=threads

  cd ${HOME}/CLEO/build/ && pwd 
  runcmd="${HOME}/CLEO/build/src/runCLEO ${configfile}"
  echo ${runcmd}
  ${runcmd}
  echo "default run"
fi 

if [[ $1 == "3" ]];
then
  export OMP_PROC_BIND=spread
  export OMP_PLACES=threads


  cd ${HOME}/CLEO/buildgpu/ && pwd 
  runcmd="${HOME}/CLEO/buildgpu/src/runCLEO ${configfile}"
  echo ${runcmd}
  ${runcmd}
  echo "gpu run"
fi