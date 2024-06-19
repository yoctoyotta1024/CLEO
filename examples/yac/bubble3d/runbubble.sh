#!/bin/bash
module load python3/2022.01-gcc-11.2.0
spack load py-numpy
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/sw/spack-levante/libfyaml-0.7.12-fvbhgo/lib
export PYTHONPATH=/work/ka1298/k202167/YAC-dev/python/:$PYTHONPATH

export OMP_PROC_BIND=spread
export OMP_PLACES=threads

mpiexec -n 1 /home/m/m300950/CLEO/build_yac/examples/yac/bubble3d/src/bubble3D \
  /home/m/m300950/CLEO/examples/yac/bubble3d/src/config/bubble3d_config.yaml \
  : -n 1 python /home/m/m300950/CLEO/examples/yac/bubble3d/yac_icon_data_reader.py
