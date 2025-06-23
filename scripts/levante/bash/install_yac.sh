#!/bin/bash
#SBATCH --job-name=install_yac
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=940M
#SBATCH --time=00:10:00
#SBATCH --account=bm1183
#SBATCH --output=./build/bin/install_yac_out.%j.out
#SBATCH --error=./build/bin/install_yac_err.%j.out

### ------------------------------------------------------- ###
### running script sucessfully installs YAC and YAXT for
### gcc 11.2.0 compiler with openmpi 4.1.2 on Levante
### and intel 2023.2.1 and openmpi 4.1.5 on Levante
### ------------------------------------------------------- ###
source /etc/profile
module purge
spack unload --all

### Note: python version used to install yac must match version used to run model
root4YAC=$1 # absolute path for YAC and YAXT installations
compilername=$2 # compile yac and yaxt with "gcc" or "intel"
python=$3 # name or absolute path to python to make YAC python bindngs with

yaxt_tag=0.11.1
yaxt_version=yaxt-${yaxt_tag}
yaxt_release_tag=release-${yaxt_tag}
yaxt_source=https://gitlab.dkrz.de/dkrz-sw/yaxt/-/archive/$yaxt_release_tag/$yaxt_version.tar.gz

yac_tag=v3.5.2
yac_version=yac_$yac_tag
yac_source=https://gitlab.dkrz.de/dkrz-sw/yac/-/archive/$yac_tag/$yac_version.tar.gz

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
bashsrc=${SCRIPT_DIR}/src
source ${bashsrc}/levante_packages.sh

if [ "${compilername}" == "" ]
then
  echo "Bad input, please specify compiler name to build yaxt and yac with"
  exit 1
elif [ "${compilername}" == "gcc" ]
then
  compiler=${levante_gcc}
  openmpi=${levante_gcc_openmpi}
  netcdf=${levante_gcc_netcdf_yac}
  netcdf_root=${levante_gcc_netcdf_root}
  fyaml_root=${levante_gcc_fyaml_root}
elif [ "${compilername}" == "intel" ]
then
  compiler=${levante_intel}
  openmpi=${levante_intel_openmpi}
  netcdf=${levante_intel_netcdf_yac}
  netcdf_root=${levante_intel_netcdf_root}
  fyaml_root=${levante_intel_fyaml_root}
else
  echo "Bad input, unrecognised compiler name"
  exit 1
fi

if [[ "${root4YAC}" == "" || "${python}" == "" ]]
then
  echo "Bad input, please specify absolute path for where you want to install YAC and python to use to make bindings"
else
  mkdir ${root4YAC}
  module load ${compiler} ${openmpi} ${netcdf}

  CC="$(command -v mpicc)"
  FC="$(command -v mpifort)"

  ### --------------------- install YAXT ------------------- ###
  mkdir ${root4YAC}/${yaxt_version}
  cd ${root4YAC}/${yaxt_version} &&  pwd
  curl -s -L ${yaxt_source} | tar xvz --strip-components=1
  mkdir build && cd build
  ../configure \
    CC=${CC} FC=${FC} \
    CFLAGS="-O0 -g -Wall" \
    FCFLAGS="-O0 -g -Wall -cpp -fimplicit-none" \
    --without-regard-for-quality \
    --without-example-programs \
    --without-perf-programs \
    --with-pic \
    --prefix=${root4YAC}/yaxt
  make -j 8
  make install
  cd ${root4YAC} && rm -rf ${yaxt_version}
  ### ------------------------------------------------------ ###

  ## --------------------- install YAC -------------------- ###
  # python bindings made in yac_version directory (note this is not yac directory!)
  mkdir ${root4YAC}/${yac_version}
  cd ${root4YAC}/${yac_version} && pwd
  curl -s -L ${yac_source} | tar xvz --strip-components=1
  mkdir build && cd build
  ../configure \
    CC=${CC} FC=${FC} \
    CFLAGS="-O0 -g -Wall" \
    FCFLAGS="-O0 -g -Wall -cpp -fimplicit-none" \
    LDFLAGS="-lm" \
    PYTHON=${python} \
    --disable-mpi-checks \
    --with-yaxt-root=${root4YAC}/yaxt \
    --with-netcdf-root=${netcdf_root} \
    --with-fyaml-root=${fyaml_root} \
    --enable-python-bindings \
    --enable-rpaths \
    --with-pic \
    --prefix=${root4YAC}/yac
  make -j 8
  make install

  mv ${root4YAC}/${yac_version}/build/python ${root4YAC}/yac/
  cd ${root4YAC} && rm -rf ${yac_version}
  ### ------------------------------------------------------ ###
fi
