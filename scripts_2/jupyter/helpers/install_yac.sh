#!/bin/bash

### ------------------------------------------------------- ###
### Running script successfully installs YAC and YAXT for
### GCC 14.3.0 with OpenMPI 5.0.8 on Jupiter.
### Note: python version used to install yac must match version used to run model.
### ------------------------------------------------------- ###

set -e
set -o pipefail # surface curl failures instead of masking them behind tar's exit status

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &>/dev/null && pwd )

root4YAC=$1                         # absolute path for YAC and YAXT installations
compilername=${2:-gcc}              # compile yac and yaxt with "gcc"
python=${3:-${CLEO_PYTHON}}         # name or absolute path to python to make YAC python bindings with

yaxt_tag=0.11.4
yaxt_version=yaxt-${yaxt_tag}
yaxt_release_tag=release-${yaxt_tag}
yaxt_source=https://gitlab.dkrz.de/dkrz-sw/yaxt/-/archive/$yaxt_release_tag/$yaxt_version.tar.gz

yac_tag=v3.9.2
yac_version=yac_$yac_tag
yac_source=https://gitlab.dkrz.de/dkrz-sw/yac/-/archive/$yac_tag/$yac_version.tar.gz

source "${SCRIPT_DIR}/jupiter_packages.sh"
jupiter_reset_modules

if [ "${compilername}" == "" ]
then
  echo "Bad input, please specify compiler name to build yaxt and yac with"
  exit 1
elif [ "${compilername}" == "gcc" ]
then
  jupiter_load_build_stack "${compilername}" "openmp"
  jupiter_load_yac_dependencies "${compilername}"
  jupiter_load_python "${compilername}" # so venv pythons linked against libpython resolve
  netcdf_root=${jupiter_gcc_netcdf_root}
  fyaml_root=${jupiter_gcc_fyaml_root}
else
  echo "Bad input, unrecognised compiler name '${compilername}'. Must be 'gcc'"
  exit 1
fi

if [[ "${root4YAC}" == "" || "${python}" == "" ]]
then
  echo "Bad input, please specify absolute path for where you want to install YAC and python to use to make bindings"
  exit 1
fi

mkdir -p ${root4YAC}
if [ ! -d ${root4YAC} ]
then
  echo "ERROR: YAC build directory not found, please make sure it exists"
  exit 1
fi

CC="$(command -v mpicc)"
FC="$(command -v mpifort)"

### --------------------- install YAXT ------------------- ###
mkdir ${root4YAC}/${yaxt_version}
cd ${root4YAC}/${yaxt_version} && pwd
curl -sS -f -L --retry 5 --retry-all-errors --retry-delay 15 ${yaxt_source} | tar xvz --strip-components=1
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
curl -sS -f -L --retry 5 --retry-all-errors --retry-delay 15 ${yac_source} | tar xvz --strip-components=1
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
