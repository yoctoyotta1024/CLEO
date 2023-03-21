#!/bin/bash

conda activate superdropsV2

CXX=/opt/homebrew/bin/g++-12 cmake -S ./ -B ./build

python ./create_gbxboundariesbinary_script.py
python ./create_initsuperdropsbinary_script.py 5e9

cd build
make clean && make

./src/coupledCVODECLEO "../src/config/config.txt" "../libs/claras_SDconstants.hpp"