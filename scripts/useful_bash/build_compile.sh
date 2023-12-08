#!/bin/bash
#SBATCH --job-name=buildCLEO
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --mem=30G
#SBATCH --time=00:05:00
#SBATCH --mail-user=clara.bayley@mpimet.mpg.de
#SBATCH --mail-type=FAIL
#SBATCH --account=mh1126
#SBATCH --output=./build/bin/buildCLEO_out.%j.out
#SBATCH --error=./build/bin/buildCLEO_err.%j.out

buildtype=$1
path2build=${HOME}/CLEO/build/

### ------------------ build_compile.sh ---------------- ###
if [ "${buildtype}" != "serial" ] && [ "${buildtype}" != "cpu" ] && [ "${buildtype}" != "gpu" ];
then
  echo "please specify the build type as 'serial', 'cpu' or 'gpu'"
fi 

if [ "${buildtype}" == "serial" ] || [ "${buildtype}" == "cpu" ] || [ "${buildtype}" == "gpu" ];
then
  echo "build type: ${buildtype}"
  echo "path to build directory: ${path2build}"

  if [[ "${buildtype}" == "serial" ]];
  then
    echo "./build_compile_serial.sh ${path2build}"
    ./build_compile_serial.sh ${path2build}

  elif [[ "${buildtype}" == "cpu" ]];
  then
    echo "./build_compile_cpus.sh ${path2build}"
    ./build_compile_cpus.sh ${path2build}

  elif [[ "${buildtype}" == "gpu" ]];
  then
    echo "./build_compile_gpus_cpus.sh ${path2build}"
    ./build_compile_gpus_cpus.sh ${path2build}
  fi
fi
### ---------------------------------------------------- ###