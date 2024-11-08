"""
Copyright (c) 2024 MPI-M, Clara Bayley


----- CLEO -----
File: bubble3d.py
Project: bubble3d
Created Date: Friday 17th November 2023
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
Script generates input files, then runs CLEO executable "bubble3d" to
piggyback ICON bubble test case
"""

import os
import shutil
import sys
from pathlib import Path

path2CLEO = Path(sys.argv[1])
path2build = Path(sys.argv[2])
config_filename = Path(sys.argv[3])
icon_grid_file = sys.argv[4]  # TODO(CB): move to config file

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
    icon_grid_file,
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
    executable = str(path2build / "examples" / "yac" / "bubble3d" / "src" / "bubble3d")
    config_filename = str(config_filename)
    print("Executable: " + str(executable))
    print("Config file: " + str(config_filename))

    cleoproc = executable + " " + config_filename
    pythonproc = path2CLEO + "/examples/yac/bubble3d/yac_icon_data_reader.py"
    cmd = "mpiexec -n 1 " + cleoproc + " : -n 1 python " + pythonproc
    os.system(cmd)


run_exectuable(path2CLEO, path2build, config_filename, dataset)
### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###

# ### ---------------------------------------------------------------- ###
# ### ------------------------- PLOT RESULTS ------------------------- ###
# ### ---------------------------------------------------------------- ###
# TODO(CB): plot results
# # read in constants and intial setup from setup .txt file
# config = pysetuptxt.get_config(setupfile, nattrs=3, isprint=True)
# consts = pysetuptxt.get_consts(setupfile, isprint=True)
# gbxs = pygbxsdat.get_gridboxes(grid_filename, consts["COORD0"], isprint=True)

# time = pyzarr.get_time(dataset)
# sddata = pyzarr.get_supers(dataset, consts)
# maxnsupers = pyzarr.get_totnsupers(dataset)
# ### ---------------------------------------------------------------- ###
# ### ---------------------------------------------------------------- ###
