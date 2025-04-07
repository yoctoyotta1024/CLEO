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
### ------------------------------------------------------- ###
source /etc/profile
module purge
spack unload --all

set -ex

root4YAC=$1 # absolute path for YAC and YAXT installations

yaxt_tag=0.11.1
yaxt_version=yaxt-${yaxt_tag}
yaxt_release_tag=release-${yaxt_tag}
yaxt_source=https://gitlab.dkrz.de/dkrz-sw/yaxt/-/archive/$yaxt_release_tag/$yaxt_version.tar.gz

yac_tag=v3.5.2
yac_version=yac_$yac_tag
yac_source=https://gitlab.dkrz.de/dkrz-sw/yac/-/archive/$yac_tag/$yac_version.tar.gz

gcc=gcc/11.2.0-gcc-11.2.0
openmpi=openmpi/4.1.2-gcc-11.2.0
netcdf=netcdf-c/4.8.1-openmpi-4.1.2-gcc-11.2.0 # must match gcc and openmpi
netcdf_root=/sw/spack-levante/netcdf-c-4.8.1-6qheqr # must match with `module show ${netcdf}`
fyaml_root=/sw/spack-levante/libfyaml-0.7.12-fvbhgo # must match with `spack location -i libfyaml``
python=python@3.9.9%gcc@=11.2.0/fwv
pycython=py-cython@0.29.33%gcc@=11.2.0/j7b4fa
pympi4py=py-mpi4py@3.1.2%gcc@=11.2.0/hdi5yl6

CC=/sw/spack-levante/openmpi-4.1.2-mnmady/bin/mpicc # must match gcc
FC=/sw/spack-levante/openmpi-4.1.2-mnmady/bin/mpif90 # must match gcc

if [ "${root4YAC}" == "" ]
then
  echo "Bad input, please specify absolute path for where you want to install YAC"
else
  mkdir ${root4YAC}
  module load ${gcc} ${openmpi} ${netcdf}
  ### ----------------- load Python ------------------------ ###
  spack load ${python}
  spack load ${pycython}
  spack load ${pympi4py}
  ### ------------------------------------------------------ ###

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
