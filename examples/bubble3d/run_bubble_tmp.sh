#!/bin/bash

### temporary file for testing running of bubble3d example
### TODO(CB): delete file, (merge with run_example like other examples)

### If you did this to your ~/.spack/upstreams.yaml to get py-netcdf from commin's:
# then python=/sw/spack-levante/python-3.9.9-fwvsvi/bin/python # is the one loaded by py-netcdf4
# upstreams:
  # community_spack:
  #   install_tree: /work/k20200/k202160/community-spack/install
  # system_installs:
  #  install_tree: /sw/spack-levante

action=$1
path2CLEO=${2:-${HOME}/CLEO}
path2yac=${3:-/work/mh1126/m300950/yac}
path2build=${4:-${HOME}/CLEO/build_bubble3d}

if [ "${action}" == "build" ]
then
  mkdir ${path2build}
  mkdir ${path2build}/bin
  mkdir ${path2build}/share

  module purge
  ${path2CLEO}/scripts/bash/build_cleo_openmp.sh \
    ${path2CLEO}/ ${path2build} true ${path2yac}

elif [ "${action}" == "compile" ]
then
  cd ${path2build}

  module purge
  module load openmpi/4.1.2-gcc-11.2.0
  module load netcdf-c/4.8.1-openmpi-4.1.2-gcc-11.2.0
  spack load openblas@0.3.18%gcc@=11.2.0

  ${path2CLEO}/scripts/bash/compile_cleo.sh \
     /work/mh1126/m300950/cleoenv openmp ${path2build} bubble3D

elif [ "${action}" == "inputfiles" ]
then
  source activate /work/mh1126/m300950/cleoenv
  /work/mh1126/m300950/cleoenv/bin/python ${path2CLEO}/examples/bubble3d/bubble3d_inputfiles.py \
   ${path2CLEO} \
   ${path2build} \
   ${path2CLEO}/examples/bubble3d/src/config/bubble3d_config.yaml \
   ${path2build}/share/bubble3d_dimlessGBxboundaries.dat \
   ${path2build}/share/bubble3d_dimlessSDsinit.dat 0

elif [ "${action}" == "run" ]
then
  cd ${path2build}

  module purge
  # note version of python must match the YAC python bindings (e.g. module load python3/2022.01-gcc-11.2.0)
  module load openmpi/4.1.2-gcc-11.2.0 # same mpi as loaded for the build
  spack load py-netcdf4

  export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/sw/spack-levante/libfyaml-0.7.12-fvbhgo/lib
  export PYTHONPATH=${PYTHONPATH}:${path2yac}/yac/python # path to python bindings

  export OMP_PROC_BIND=spread
  export OMP_PLACES=threads

  cp /work/mh1126/m300950/icon/build/experiments/aes_bubble/aes_bubble_atm_3d_ml_20080801T000000Z.nc \
    ${path2build}/aes_bubble_atm_3d_ml_20080801T000000Z.nc

  cp /work/mh1126/m300950/icon/build/experiments/aes_bubble/Torus_Triangles_20x4_5000m.nc \
    ${path2build}/Torus_Triangles_20x4_5000m.nc

  mpiexec -n 1 ${path2build}/examples/bubble3d/src/bubble3D \
    ${path2CLEO}/examples/bubble3d/src/config/bubble3d_config.yaml \
    : -n 1 python \
    ${path2CLEO}/examples/bubble3d/yac_bubble_data_reader.py \
    ${path2build}/aes_bubble_atm_3d_ml_20080801T000000Z.nc \
    ${path2build}/Torus_Triangles_20x4_5000m.nc
fi
