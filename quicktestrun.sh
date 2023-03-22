#!/bin/bash

module load python3/2022.01-gcc-11.2.0
source activate /work/mh1126/m300950/pySDenv

python ./create_gbxboundariesbinary_script.py ./
python ./create_initsuperdropsbinary_script.py ./ 1e9

cd build
make clean && make
./src/coupledCVODECLEO "../src/config/config.txt" "../libs/claras_SDconstants.hpp"