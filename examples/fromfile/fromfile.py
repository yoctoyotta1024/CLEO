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
import argparse
from pathlib import Path
import fromfile_inputfiles
import fromfile_plotting

parser = argparse.ArgumentParser()
parser.add_argument(
    "path2CLEO", type=Path, help="Absolute path to CLEO directory (for PySD)"
)
parser.add_argument("path2build", type=Path, help="Absolute path to build directory")
parser.add_argument(
    "config_filename", type=Path, help="Absolute path to configuration YAML file"
)
parser.add_argument(
    "--ntasks",
    type=int,
    default="4",
    help="Number of MPI processes to run program with",
)
parser.add_argument(
    "--do_inputfiles",
    type=str,
    choices=["TRUE", "FALSE"],
    default="TRUE",
    help="Generate initial condition binary files",
)
parser.add_argument(
    "--do_run_executable",
    type=str,
    choices=["TRUE", "FALSE"],
    default="TRUE",
    help="Run fromfile executable",
)
parser.add_argument(
    "--do_plot_results",
    type=str,
    choices=["TRUE", "FALSE"],
    default="TRUE",
    help="Plot results of fromfile example",
)
args = parser.parse_args()

path2CLEO = args.path2CLEO
path2build = args.path2build
config_filename = args.config_filename
ntasks = args.ntasks

do_inputfiles = True
if args.do_inputfiles == "FALSE":
    do_inputfiles = False
do_run_executable = True
if args.do_run_executable == "FALSE":
    do_run_executable = False
do_plot_results = True
if args.do_plot_results == "FALSE":
    do_plot_results = False

isfigures = [True, True]  # booleans for [making, saving] initialisation figures

### ---------------------------------------------------------------- ###
### ----------------------- INPUT PARAMETERS ----------------------- ###
### ---------------------------------------------------------------- ###
### --- essential paths and filenames --- ###
# path and filenames for creating initial SD conditions
binpath = path2build / "bin" / f"ntasks{ntasks}"
sharepath = path2build / "share"
grid_filename = sharepath / "fromfile_dimlessGBxboundaries.dat"
initsupers_filename = sharepath / "fromfile_dimlessSDsinit.dat"
thermofiles = sharepath / "fromfile_dimlessthermo.dat"
savefigpath = binpath  # directory for saving figures

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
        binpath.parent.mkdir(exist_ok=True)
        binpath.mkdir(exist_ok=True)
        savefigpath.parent.mkdir(exist_ok=True)
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
        savefigpath,
        config_filename,
        grid_filename,
        initsupers_filename,
        thermofiles,
        isfigures=isfigures,
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
    subprocess.run(["srun", f"--ntasks={ntasks}", executable, config_filename])
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
