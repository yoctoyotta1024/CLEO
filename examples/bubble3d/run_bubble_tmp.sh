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
path2yac=${3:-/work/bm1183/m300950/yacyaxt}
path2build=${4:-${HOME}/CLEO/build_bubble3d}

icon_grid_file=/work/bm1183/m300950/icon/build/experiments/aes_bubble/aes_bubble_atm_cgrid_ml.nc
icon_data_file=/work/bm1183/m300950/icon/build/experiments/aes_bubble/aes_bubble_atm_3d_ml_20080801T000000Z.nc
icon_grid_file_copy=${path2build}/share/icon_grid_file_aes_bubble_atm_cgrid_ml.nc
icon_data_file_copy=${path2build}/share/icon_data_file_aes_bubble_atm_3d_ml_20080801T000000Z.nc

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
  spack unload --all
  module load openmpi/4.1.2-gcc-11.2.0
  module load netcdf-c/4.8.1-openmpi-4.1.2-gcc-11.2.0
  spack load openblas@0.3.18%gcc@=11.2.0

  ${path2CLEO}/scripts/bash/compile_cleo.sh openmp ${path2build} bubble3d

elif [ "${action}" == "inputfiles" ]
then
  cp ${icon_grid_file} ${icon_grid_file_copy}
  cp ${icon_data_file} ${icon_data_file_copy}

  cleoenv=/work/bm1183/m300950/bin/envs/cleoenv
  python=${cleoenv}/bin/python3
  ${python} ${path2CLEO}/examples/bubble3d/bubble3d_inputfiles.py \
   ${path2CLEO} \
   ${path2build} \
   ${path2CLEO}/examples/bubble3d/src/config/bubble3d_config.yaml \
   ${path2build}/share/bubble3d_dimlessGBxboundaries.dat \
   ${path2build}/share/bubble3d_dimlessSDsinit.dat \
   ${icon_grid_file_copy} \
   0

elif [ "${action}" == "run" ]
then
  cd ${path2build}

  module purge
  spack unload --all
  # note version of python must match the YAC python bindings (e.g. module load python3/2022.01-gcc-11.2.0)
  module load openmpi/4.1.2-gcc-11.2.0 # same mpi as loaded for the build
  ### ----------------- load Python ------------------------ ###
  spack load python@3.9.9%gcc@=11.2.0/fwv
  spack load py-cython@0.29.33%gcc@=11.2.0/j7b4fa
  spack load py-mpi4py@3.1.2%gcc@=11.2.0/hdi5yl6
  ### ------------------------------------------------------ ###

  export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/sw/spack-levante/libfyaml-0.7.12-fvbhgo/lib
  export PYTHONPATH=${PYTHONPATH}:${path2yac}/yac/python # path to python bindings

  export OMP_PROC_BIND=spread
  export OMP_PLACES=threads

  icon_grid_name="icon_atmos_grid" # must match CLEO (see yac_cartesian_dynamics.cpp)
  icon_data_timestep=30 # must match ICON data file [seconds]
  cleo_coupling_timestep=60 # must match CLEO config file [seconds]
  cleo_t_end=7200 # must match CLEO config file [seconds]
  cleo_num_vertical_levels=24 # must match CLEO grid_filename

  mpiexec -n 1 ${path2build}/examples/bubble3d/src/bubble3d \
    ${path2CLEO}/examples/bubble3d/src/config/bubble3d_config.yaml \
    : -n 1 python \
    ${path2CLEO}/examples/bubble3d/yac_bubble_data_reader.py \
    ${icon_grid_file_copy} ${icon_data_file_copy} ${icon_grid_name} ${icon_data_timestep} \
    ${cleo_coupling_timestep} ${cleo_t_end} ${cleo_num_vertical_levels}
fi
