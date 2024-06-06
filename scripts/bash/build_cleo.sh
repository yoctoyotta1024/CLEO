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

buildtype=$1  # required
path2CLEO=$2  # required
path2build=$3 # required
enableyac=$4  # required "true" or otherwise
yacyaxtroot=$5    # required if enableyac == "true"

path2buildbash=${path2CLEO}/scripts/bash/

if [[ "${buildtype}" == "" ||
      "${path2CLEO}" == "" ||
      "${path2build}" == "" ||
      "${path2CLEO}" == "${path2build}" ]]
then
  echo "Bad inputs, please check your buildtype, path2CLEO and path2build"
  exit 0

else
  # check that path to installation is set if YAC is enabled
  if [[ ${enableyac} == "true" && "${yacyaxtroot}"  == "" ]]
    then
      echo "Bad YAC flags, if enableyac==true please provide root to YAC installation"
      exit 0
  else
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

      if [ ${enableyac} == "true" ]
      then
        echo "Build with YAC enabled from ${yacyaxtroot}/yac and ${yacyaxtroot}/yaxt"
      else
        echo "Build with YAC disabled"
      fi

      if [[ "${buildtype}" == "serial" ]];
      then
        echo "${path2buildbash}/build_cleo_serial.sh ${path2CLEO} ${path2build}"
        ${path2buildbash}/build_cleo_serial.sh ${path2CLEO} ${path2build} ${enableyac} ${yacyaxtroot}

      elif [[ "${buildtype}" == "openmp" ]];
      then
      echo "${path2buildbash}/build_cleo_openmp.sh ${path2CLEO} ${path2build}"
        ${path2buildbash}/build_cleo_openmp.sh ${path2CLEO} ${path2build} ${enableyac} ${yacyaxtroot}

      elif [[ "${buildtype}" == "cuda" ]];
      then
        echo "${path2buildbash}/build_cleo_cuda_openmp.sh ${path2CLEO} ${path2build}"
        ${path2buildbash}/build_cleo_cuda_openmp.sh ${path2CLEO} ${path2build} ${enableyac} ${yacyaxtroot}
      fi
    fi
    ### ---------------------------------------------------- ###
  fi
fi
