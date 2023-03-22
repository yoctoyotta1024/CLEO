#!/bin/bash

module load gcc/11.2.0-gcc-11.2.0
conda activate /work/mh1126/m300950/superdropsV2

CXX=g++ cmake -S ./ -B ./build

python ./create_gbxboundariesbinary_script.py ./
python ./create_initsuperdropsbinary_script.py ./ 5e9

cd build
make clean && make

./src/coupledCVODECLEO "../src/config/config.txt" "../libs/claras_SDconstants.hpp"