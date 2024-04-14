#!/bin/bash
#SBATCH --job-name=shima2009
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --mem=30G
#SBATCH --time=00:05:00
#SBATCH --mail-user=clara.bayley@mpimet.mpg.de
#SBATCH --mail-type=FAIL
#SBATCH --account=mh1126
#SBATCH --output=./shima2009_out.%j.out
#SBATCH --error=./shima2009_err.%j.out

### ---------------------------------------------------- ###
### ------------------ Input Parameters ---------------- ###
### ------ You MUST edit these lines to set your ------- ###
### ---- build type, directories, the executable(s) ---- ###
### -------- to compile, and your python script -------- ###
### ---------------------------------------------------- ###
buildtype="cuda"
path2CLEO=${HOME}/CLEO/
path2build=${HOME}/CLEO/build_shima2009/
executables="golcolls longcolls lowlistcolls"

pythonscript=${path2CLEO}/examples/boxmodelcollisions/shima2009.py
configfile=${path2CLEO}/examples/boxmodelcollisions/shima2009_config.txt
script_args=${configfile}" golovin long lowlist"
### ---------------------------------------------------- ###
### ---------------------------------------------------- ###
### ---------------------------------------------------- ###

### ---------- build, compile and run example ---------- ###
runcmd="${path2CLEO}/examples/run_example.sh ${buildtype} \
  ${path2CLEO} ${path2build} ${executables} ${pythonscript} ${script_args}"
echo ${runcmd}
${runcmd}
### ---------------------------------------------------- ###
