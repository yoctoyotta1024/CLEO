"""
Copyright (c) 2024 MPI-M, Clara Bayley


----- CLEO -----
File: bubble3d.py
Project: bubble3d
Created Date: Friday 17th November 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
Last Modified: Tuesday 15th April 2025
Modified By: CB
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
Script generates input files, then runs CLEO executable "bubble3d" to
piggyback ICON bubble test case
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path

path2CLEO = Path(sys.argv[1])
path2build = Path(sys.argv[2])
config_filename = Path(sys.argv[3])

import bubble3d_inputfiles

### ---------------------------------------------------------------- ###
### ----------------------- INPUT PARAMETERS ----------------------- ###
### ---------------------------------------------------------------- ###
### --- essential paths and filenames --- ###
# path and filenames for creating initial SD conditions
binpath = path2build / "bin"
sharepath = path2build / "share"
grid_filename = sharepath / "bubble3d_dimlessGBxboundaries.dat"
initsupers_filename = sharepath / "bubble3d_dimlessSDsinit.dat"
savefigpath = path2build / "bin"
SDgbxs2plt = [0]  # gbxindex of SDs to plot (nb. "all" can be very slow)

# path and file names for plotting results
setupfile = binpath / "bubble3d_setup.txt"
dataset = binpath / "bubble3d_sol.zarr"

### ---------------------------------------------------------------- ###
### ------------------- BINARY FILES GENERATION--------------------- ###
### ---------------------------------------------------------------- ###
### --- ensure build, share and bin directories exist --- ###
if path2CLEO == path2build:
    raise ValueError("build directory cannot be CLEO")
else:
    path2build.mkdir(exist_ok=True)
    sharepath.mkdir(exist_ok=True)
    binpath.mkdir(exist_ok=True)
    savefigpath.mkdir(exist_ok=True)

### --- delete any existing initial conditions --- ###
shutil.rmtree(grid_filename, ignore_errors=True)
shutil.rmtree(initsupers_filename, ignore_errors=True)

bubble3d_inputfiles.main(
    path2CLEO,
    path2build,
    config_filename,
    grid_filename,
    initsupers_filename,
    SDgbxs2plt,
)
### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###


### ---------------------------------------------------------------- ###
### ---------------------- RUN CLEO EXECUTABLE --------------------- ###
### ---------------------------------------------------------------- ###
def run_exectuable(path2CLEO, path2build, config_filename, dataset):
    """delete existing dataset, the run exectuable with given config file"""
    os.chdir(path2build)
    os.system("pwd")
    shutil.rmtree(dataset, ignore_errors=True)  # delete any existing dataset
    cleoproc = str(path2build / "examples" / "bubble3d" / "src" / "bubble3d")
    cleoproc_args = str(config_filename)
    print("CLEO Executable: " + cleoproc)
    print("CLEO Config file: " + cleoproc_args)

    python = sys.executable
    pythonproc = str(path2CLEO / "examples" / "bubble3d" / "yac_bubble_data_reader.py")
    pythonproc_args = [str(path2build), str(config_filename)]
    print("YAC script: " + pythonproc)
    print("YAC arguments: " + " ".join(pythonproc_args))

    cmd = [
        "mpiexec",
        "-n",
        "1",
        cleoproc,
        cleoproc_args,
        ":",
        "-n",
        "1",
        python,
        pythonproc,
    ] + pythonproc_args
    print(" ".join(cmd))
    subprocess.run(cmd)


run_exectuable(path2CLEO, path2build, config_filename, dataset)
### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###


### ---------------------------------------------------------------- ###
### ------------------------- PLOT RESULTS ------------------------- ###
### ---------------------------------------------------------------- ###
def plot_results(path2CLEO):
    plotting_script = path2CLEO / "examples" / "bubble3d" / "bubble3d_plotting.py"
    python = sys.executable
    cmd = [python, plotting_script, path2CLEO, path2build, config_filename]
    subprocess.run(cmd)


plot_results(path2CLEO)
### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###
