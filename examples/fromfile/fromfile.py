"""
Copyright (c) 2024 MPI-M, Clara Bayley


----- CLEO -----
File: fromfile.py
Project: fromfile
Created Date: Friday 17th November 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
Last Modified: Tuesday 9th July 2024
Modified By: CB
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
Script generates input files, then runs CLEO executable for 3D example with
time varying thermodynamics read from binary files. Plots output data as .png
files for visual checks.
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path
import fromfile_inputfiles
import fromfile_plotting

path2CLEO = Path(sys.argv[1])
path2build = Path(sys.argv[2])
config_filename = Path(sys.argv[3])

### ---------------------------------------------------------------- ###
### ----------------------- INPUT PARAMETERS ----------------------- ###
### ---------------------------------------------------------------- ###
### --- essential paths and filenames --- ###
do_inputfiles = True
do_run_executable = True
do_plot_results = True

### --- essential paths and filenames --- ###
# path and filenames for creating initial SD conditions
binpath = path2build / "bin"
sharepath = path2build / "share"
grid_filename = sharepath / "fromfile_dimlessGBxboundaries.dat"
initsupers_filename = sharepath / "fromfile_dimlessSDsinit.dat"
thermofiles = sharepath / "fromfile_dimlessthermo.dat"
savefigpath = path2build / "bin"  # directory for saving figures

# path and file names for plotting results
setupfile = binpath / "fromfile_setup.txt"
dataset = binpath / "fromfile_sol.zarr"
### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###

### ---------------------------------------------------------------- ###
### ------------------- BINARY FILES GENERATION--------------------- ###
### ---------------------------------------------------------------- ###
if do_inputfiles:
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
    shutil.rmtree(
        str(thermofiles)[:-4] + "*", ignore_errors=True
    )  # delete any existing dataset

    fromfile_inputfiles.main(
        path2CLEO,
        path2build,
        config_filename,
        grid_filename,
        initsupers_filename,
        thermofiles,
    )
### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###

### ---------------------------------------------------------------- ###
### ---------------------- RUN CLEO EXECUTABLE --------------------- ###
### ---------------------------------------------------------------- ###
if do_run_executable:
    os.chdir(path2build)
    subprocess.run(["pwd"])
    shutil.rmtree(dataset, ignore_errors=True)  # delete any existing dataset
    executable = path2build / "examples" / "fromfile" / "src" / "fromfile"
    print("Executable: " + str(executable))
    print("Config file: " + str(config_filename))
    subprocess.run(["srun", "--ntasks=4", executable, config_filename])
### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###

### ---------------------------------------------------------------- ###
### ------------------------- PLOT RESULTS ------------------------- ###
### ---------------------------------------------------------------- ###
if do_plot_results:
    fromfile_plotting.main(
        path2CLEO,
        grid_filename,
        setupfile,
        dataset,
        savefigpath,
    )
### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###
