#!/bin/bash

### First add this to your ~/.spack/upstreams.yaml to get py-netcdf from commin's:
# upstreams:
  # community_spack:
  #   install_tree: /work/k20200/k202160/community-spack/install
  # system_installs:
  #  install_tree: /sw/spack-levante

# Note: python=/sw/spack-levante/python-3.9.9-fwvsvi/bin/python # the one loaded by py-netcdf4

cd /home/m/m300950/CLEO/build_yac

# spack load py-netcdf4
module load python3/2022.01-gcc-11.2.0
source activate /work/mh1126/m300950/cleoenv
module load openmpi/4.1.2-gcc-11.2.0 # same mpi as loaded for the build

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/sw/spack-levante/libfyaml-0.7.12-fvbhgo/lib
export PYTHONPATH=/work/mh1126/m300950/yac/yac/python/:$PYTHONPATH

export OMP_PROC_BIND=spread
export OMP_PLACES=threads

/work/mh1126/m300950/cleoenv/bin/mpiexec -n 1 /home/m/m300950/CLEO/build_yac/examples/yac/bubble3d/src/bubble3D \
  /home/m/m300950/CLEO/examples/yac/bubble3d/src/config/bubble3d_config.yaml \
  : -n 1 /work/mh1126/m300950/cleoenv/bin/python3 \
  /home/m/m300950/CLEO/examples/yac/bubble3d/yac_icon_data_reader.py
