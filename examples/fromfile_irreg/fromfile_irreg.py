"""
Copyright (c) 2024 MPI-M, Clara Bayley


----- CLEO -----
File: fromfile_irreg.py
Project: fromfile_irreg
Created Date: Wednesday 11th September 2024
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
Script generates input files, then runs CLEO executable for 3D example with
irreglar 3D grid and time varying thermodynamics read from binary files.
Plots output data as .png files for visual checks.
"""

import os
import shutil
import subprocess
import argparse
from pathlib import Path
import fromfile_irreg_inputfiles
import fromfile_irreg_plotting

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
    action="store_true",  # default is False
    help="Generate initial condition binary files",
)
parser.add_argument(
    "--do_run_executable",
    action="store_true",  # default is False
    help="Run fromfile executable",
)
parser.add_argument(
    "--do_plot_results",
    action="store_true",  # default is False
    help="Plot results of fromfile example",
)
args = parser.parse_args()

path2CLEO = args.path2CLEO
path2build = args.path2build
config_filename = args.config_filename
ntasks = args.ntasks

isfigures = [False, True]  # booleans for [showing, saving] initialisation figures

### ---------------------------------------------------------------- ###
### ----------------------- INPUT PARAMETERS ----------------------- ###
### ---------------------------------------------------------------- ###
### --- essential paths and filenames --- ###
# path and filenames for creating initial SD conditions
binpath = path2build / "bin" / f"ntasks{ntasks}"
sharepath = path2build / "share"
grid_filename = sharepath / "fromfile_irreg_dimlessGBxboundaries.dat"
initsupers_filename = sharepath / "fromfile_irreg_dimlessSDsinit.dat"
thermofiles = sharepath / "fromfile_irreg_dimlessthermo.dat"
savefigpath = binpath  # directory for saving figures

# path and file names for plotting results
setupfile = binpath / "fromfile_irreg_setup.txt"
dataset = binpath / "fromfile_irreg_sol.zarr"
### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###

### ---------------------------------------------------------------- ###
### ------------------- BINARY FILES GENERATION--------------------- ###
### ---------------------------------------------------------------- ###
if args.do_inputfiles:
    ### --- ensure build, share and bin directories exist --- ###
    if path2CLEO == path2build:
        raise ValueError("build directory cannot be CLEO")
    else:
        path2build.mkdir(exist_ok=True)
        sharepath.mkdir(exist_ok=True)
        binpath.parent.mkdir(exist_ok=True)
        binpath.mkdir(exist_ok=True)
        savefigpath.mkdir(exist_ok=True)

    ### --- delete any existing initial conditions --- ###
    shutil.rmtree(grid_filename, ignore_errors=True)
    shutil.rmtree(initsupers_filename, ignore_errors=True)
    all_thermofiles = thermofiles.parent / Path(
        f"{thermofiles.stem}*{thermofiles.suffix}"
    )
    shutil.rmtree(all_thermofiles, ignore_errors=True)

    fromfile_irreg_inputfiles.main(
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
if args.do_run_executable:
    os.chdir(path2build)
    subprocess.run(["pwd"])
    shutil.rmtree(dataset, ignore_errors=True)  # delete any existing dataset
    executable = path2build / "examples" / "fromfile_irreg" / "src" / "fromfile_irreg"
    print("Executable: " + str(executable))
    print("Config file: " + str(config_filename))
    subprocess.run(["srun", f"--ntasks={ntasks}", executable, config_filename])
### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###

### ---------------------------------------------------------------- ###
### ------------------------- PLOT RESULTS ------------------------- ###
### ---------------------------------------------------------------- ###
if args.do_plot_results:
    fromfile_irreg_plotting.main(
        path2CLEO,
        grid_filename,
        setupfile,
        dataset,
        savefigpath,
    )
### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###
