#!/bin/bash
#SBATCH --job-name=build_cleo
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --mem=30G
#SBATCH --time=00:05:00
#SBATCH --mail-user=clara.bayley@mpimet.mpg.de
#SBATCH --mail-type=FAIL
#SBATCH --account=mh1126
#SBATCH --output=./build/bin/build_cleo_out.%j.out
#SBATCH --error=./build/bin/build_cleo_err.%j.out

buildtype=$1
path2CLEO=$2
path2build=$3
path2buildbash=${path2CLEO}/scripts/bash/

### --------------------- build CLEO ------------------- ###
if [ "${buildtype}" != "serial" ] && [ "${buildtype}" != "openmp" ] && [ "${buildtype}" != "cuda" ];
then
  echo "please specify the build type as 'serial', 'openmp' or 'cuda'"
fi

if [ "${buildtype}" == "serial" ] || [ "${buildtype}" == "openmp" ] || [ "${buildtype}" == "cuda" ];
then
  echo "build type: ${buildtype}"
  echo "path to CLEO: ${path2CLEO}"
  echo "path to build directory: ${path2build}"

  if [[ "${buildtype}" == "serial" ]];
  then
    echo "${path2buildbash}/build_cleo_serial.sh ${path2CLEO} ${path2build}"
    ${path2buildbash}/build_cleo_serial.sh ${path2CLEO} ${path2build}

  elif [[ "${buildtype}" == "openmp" ]];
  then
   echo "${path2buildbash}/build_cleo_openmp.sh ${path2CLEO} ${path2build}"
    ${path2buildbash}/build_cleo_openmp.sh ${path2CLEO} ${path2build}

  elif [[ "${buildtype}" == "cda" ]];
  then
    echo "${path2buildbash}/build_cleo_openmp_cuda.sh ${path2CLEO} ${path2build}"
    ${path2buildbash}/build_cleo_openmp_cuda.sh ${path2CLEO} ${path2build}
  fi
fi
### ---------------------------------------------------- ###
