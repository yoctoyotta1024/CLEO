#!/bin/bash

### temporary file for testing running of bubble3d example
### TODO(CB): delete file

### If you did this to your ~/.spack/upstreams.yaml to get py-netcdf from commin's:
# then python=/sw/spack-levante/python-3.9.9-fwvsvi/bin/python # is the one loaded by py-netcdf4
# upstreams:
  # community_spack:
  #   install_tree: /work/k20200/k202160/community-spack/install
  # system_installs:
  #  install_tree: /sw/spack-levante

source /etc/profile

action=$1
path2CLEO=${2:-${HOME}/CLEO}
path2yac=${3:-/work/bm1183/m300950/yacyaxt/gcc}
path2build=${4:-${HOME}/CLEO/build_bubble3d}
python=${5:-/home/m/m300950/CLEO/.venv/bin/python3}

if [ "${action}" == "build" ]
then
  mkdir ${path2build}
  mkdir ${path2build}/bin
  mkdir ${path2build}/share

  gcc=gcc/11.2.0-gcc-11.2.0
  openmpi=openmpi/4.1.2-gcc-11.2.0
  netcdf=netcdf-c/4.8.1-openmpi-4.1.2-gcc-11.2.0 # must match gcc and openmpi
  spack load openblas@0.3.18%gcc@=11.2.0
  module load ${gcc} ${openmpi} ${netcdf}

  cmake -S ${path2CLEO} -B ${path2build} -DCLEO_COUPLED_DYNAMICS="yac" \
    -DCLEO_YAC_MODULE_PATH="${path2CLEO}/libs/coupldyn_yac/cmake" \
    -DCLEO_YAXT_ROOT=${path2yac}/yaxt -DCLEO_YAC_ROOT=${path2yac}/yac

elif [ "${action}" == "compile" ]
then
  module purge
  spack unload --all
  module load openmpi/4.1.2-gcc-11.2.0

  cd ${path2build}
  make -j 128 bubble3d

elif [ "${action}" == "inputfiles" ]
then
  ${python} ${path2CLEO}/examples/bubble3d/bubble3d_inputfiles.py \
   ${path2CLEO} \
   ${path2build} \
   ${path2CLEO}/examples/bubble3d/src/config/bubble3d_config.yaml \
   --save_figures \
   --savefigpath="${path2build}/bin"

elif [ "${action}" == "run" ]
then
  module purge
  spack unload --all
  # note version of python must match the YAC python bindings (e.g. spack load python@3.9.9%gcc@=11.2.0/fwv)
  module load openmpi/4.1.2-gcc-11.2.0 # same mpi as loaded for the build

  export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/sw/spack-levante/libfyaml-0.7.12-fvbhgo/lib
  export PYTHONPATH=${PYTHONPATH}:${path2yac}/yac/python # path to python bindings

  export OMP_PROC_BIND=spread
  export OMP_PLACES=threads

  cd ${path2build}/../
  mpiexec -n 1 ${path2build}/examples/bubble3d/src/bubble3d \
    ${path2CLEO}/examples/bubble3d/src/config/bubble3d_config.yaml \
    : -n 1 ${python} \
    ${path2CLEO}/examples/bubble3d/yac_bubble_data_reader.py \
    ${path2build} \
    ${path2CLEO}/examples/bubble3d/src/config/bubble3d_config.yaml

elif [ "${action}" == "plot" ]
then
  ${python} ${path2CLEO}/examples/bubble3d/bubble3d_plotting.py \
    ${path2CLEO} \
    ${path2build}/bin \
    ${path2build}/share/bubble3d_dimlessGBxboundaries.dat \
    ${path2build}/bin/bubble3d_setup.txt \
    ${path2build}/bin/bubble3d_sol.zarr
fi
