"""
Copyright (c) 2024 MPI-M, Clara Bayley


----- CLEO -----
File: fromfile_irreg.py
Project: fromfile_irreg
Created Date: Wednesday 11th September 2024
Author: Clara Bayley (CB)
Additional Contributors:
-----
Last Modified: Wednesday 11th September 2024
Modified By: CB
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
import sys
from pathlib import Path
import fromfile_irreg_inputfiles

path2CLEO = Path(sys.argv[1])
path2build = Path(sys.argv[2])
config_filename = Path(sys.argv[3])

sys.path.append(str(path2CLEO))  # imports from pySD
sys.path.append(
    str(path2CLEO / "examples" / "exampleplotting")
)  # imports from example plots package

from src import plot_output_thermo
from plotssrc import pltsds, pltmoms
from pySD.sdmout_src import pyzarr, pysetuptxt, pygbxsdat

### ---------------------------------------------------------------- ###
### ----------------------- INPUT PARAMETERS ----------------------- ###
### ---------------------------------------------------------------- ###
### --- essential paths and filenames --- ###
# path and filenames for creating initial SD conditions
binpath = path2build / "bin"
sharepath = path2build / "share"
grid_filename = sharepath / "fromfile_irreg_dimlessGBxboundaries.dat"
initsupers_filename = sharepath / "fromfile_irreg_dimlessSDsinit.dat"
thermofiles = sharepath / "fromfile_irreg_dimlessthermo.dat"
savefigpath = path2build / "bin"  # directory for saving figures

# path and file names for plotting results
setupfile = binpath / "fromfile_irreg_setup.txt"
dataset = binpath / "fromfile_irreg_sol.zarr"
### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###

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
all_thermofiles = thermofiles.parent / Path(f"{thermofiles.stem}*{thermofiles.suffix}")
shutil.rmtree(all_thermofiles, ignore_errors=True)

fromfile_irreg_inputfiles.main(
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
os.chdir(path2build)
subprocess.run(["pwd"])
shutil.rmtree(dataset, ignore_errors=True)  # delete any existing dataset
executable = path2build / "examples" / "fromfile_irreg" / "src" / "fromfile_irreg"
print("Executable: " + str(executable))
print("Config file: " + str(config_filename))
subprocess.run([executable, config_filename])
### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###


### ---------------------------------------------------------------- ###
### ------------------------- PLOT RESULTS ------------------------- ###
### ---------------------------------------------------------------- ###
# read in constants and intial setup from setup .txt file
config = pysetuptxt.get_config(setupfile, nattrs=3, isprint=True)
consts = pysetuptxt.get_consts(setupfile, isprint=True)
gbxs = pygbxsdat.get_gridboxes(grid_filename, consts["COORD0"], isprint=True)

time = pyzarr.get_time(dataset)
sddata = pyzarr.get_supers(dataset, consts)
maxnsupers = pyzarr.get_totnsupers(dataset)
thermo, winds = pyzarr.get_thermodata(
    dataset, config["ntime"], gbxs["ndims"], consts, getwinds=True
)

# plot super-droplet results
savename = savefigpath / "fromfile_irreg_maxnsupers_validation.png"
pltmoms.plot_totnsupers(time, maxnsupers, savename=savename)

nsample = 1000
savename = savefigpath / "fromfile_irreg_motion2d_validation.png"
pltsds.plot_randomsample_superdrops_2dmotion(
    sddata,
    config["maxnsupers"],
    nsample,
    savename=savename,
    arrows=False,
    israndom=False,
)

# plot thermodynamics results
plot_output_thermo.plot_domain_thermodynamics_timeseries(
    time, gbxs, thermo, winds, savedir=savefigpath
)

### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###
