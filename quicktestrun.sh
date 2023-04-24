#!/bin/bash

module load python3/2022.01-gcc-11.2.0
source activate /work/mh1126/m300950/pySDenv

# cd build
# make clean && make
./src/runCLEO "../src/config/config.txt" "../libs/claras_SDconstants.hpp"