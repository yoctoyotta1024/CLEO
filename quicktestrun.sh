#!/bin/bash
#SBATCH --job-name=quickrun
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --mem=30G
#SBATCH --time=00:30:00
#SBATCH --mail-user=clara.bayley@mpimet.mpg.de
#SBATCH --mail-type=FAIL
#SBATCH --account=mh1126
#SBATCH --output=./build/bin/quickrun_out.%j.out
#SBATCH --error=./build/bin/quickrun_err.%j.out

module load python3/2022.01-gcc-11.2.0
source activate /work/mh1126/m300950/pySDenv
module load gcc/11.2.0-gcc-11.2.0

cd build
# make clean && make
#./src/runCLEO "../src/config/config.txt" "../libs/claras_SDconstants.hpp"
./src/cond0D "../src/config/condconfig.txt" "../libs/claras_SDconstants.hpp"