#!/bin/bash
#SBATCH --job-name=installyac
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --mem=30G
#SBATCH --time=00:10:00
#SBATCH --mail-user=clara.bayley@mpimet.mpg.de
#SBATCH --mail-type=FAIL
#SBATCH --account=mh1126
#SBATCH --output=./build/bin/installyac_out.%j.out
#SBATCH --error=./build/bin/installyac_err.%j.out

### ------------------------------------------------------- ###
### running script sucessfully installs YAC and YAXT for
### gcc 11.2.0 compiler with openmpi 4.1.2 on Levante
### ------------------------------------------------------- ###

root4YAC=$1 # absolute path for YAC and YAXT installations

yaxt_source=https://swprojects.dkrz.de/redmine/attachments/download/534/yaxt-0.11.1.tar.gz
yaxt_version=yaxt-0.11.1 # must match yaxt_source

yac_source=https://gitlab.dkrz.de/dkrz-sw/yac/-/archive/v3.2.0_b_p1/yac-v3.2.0_b_p1.tar.gz
yac_version=yac-v3.2.0_b_p1 # must match yac_source

gcc=gcc/11.2.0-gcc-11.2.0
openmpi=openmpi/4.1.2-gcc-11.2.0
hdf5=hdf5/1.12.1-openmpi-4.1.2-gcc-11.2.0
inteloneapi=intel-oneapi-mkl/2022.0.1-gcc-11.2.0
netcdf=netcdf-c/4.8.1-openmpi-4.1.2-gcc-11.2.0 # must match gcc and openmpi
netcdf_root=/sw/spack-levante/netcdf-c-4.8.1-6qheqr # must match with `module show ${netcdf}`
fyaml=libfyaml
fyaml_root=/sw/spack-levante/libfyaml-0.7.12-fvbhgo # must match with `spack location -i libfyaml``

CC=/sw/spack-levante/openmpi-4.1.2-mnmady/bin/mpicc # must match gcc
FC=/sw/spack-levante/openmpi-4.1.2-mnmady/bin/mpif90 # must match gcc

if [ "${root4YAC}" == "" ]
then
  echo "Bad input, please specify absolute path for where you want to install YAC"
else
  module load ${gcc} ${openmpi} ${netcdf} ${hdf5} ${inteloneapi}
  spack load ${fyaml}

  ### --------------------- install YAXT ------------------- ###
  cd ${root4YAC} && pwd
  curl -s -L ${yaxt_source} | tar xvz
  cd ${yaxt_version}
  ./configure \
    CC=${CC} FC=${FC} \
    CFLAGS="-O0 -g -Wall -fPIC" \
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
  cd ${root4YAC} && pwd
  curl -s -L ${yac_source} | tar xvz
  cd ${yac_version}
  ./configure \
    CC=${CC} FC=${FC} \
    CFLAGS="-O0 -g -Wall -fPIC" \
    FCFLAGS="-O0 -g -Wall -cpp -fimplicit-none" \
    MKL_CLIBS="`pkg-config --libs mkl-static-lp64-seq`" \
    MKL_CFLAGS="`pkg-config --cflags mkl-static-lp64-seq`" \
    LDFLAGS="-lm" \
    --disable-mpi-checks \
    --with-yaxt-root=${root4YAC}/yaxt \
    --with-netcdf-root=${netcdf_root} \
    --with-fyaml-root=${fyaml_root} \
    --enable-rpaths \
    --prefix=${root4YAC}/yac
  make -j 8
  make install
  cd ${root4YAC} && rm -rf ${yac_version}
### ------------------------------------------------------ ###
fi
