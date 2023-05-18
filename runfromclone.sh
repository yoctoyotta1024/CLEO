#!/bin/bash

module load gcc/11.2.0-gcc-11.2.0

module load python3/2022.01-gcc-11.2.0
source activate /work/mh1126/m300950/pySDenv

CXX=g++ cmake -S ./ -B ./build

python ./create_gbxboundariesbinary_script.py ./
python ./create_initsuperdropsbinary_script.py ./ 
python ./create_initthermobinary_script.py ./

cd build
make clean && make -j 16
./src/runCLEO "../src/config/config.txt" "../libs/claras_SDconstants.hpp"