"""
Copyright (c) 2024 MPI-M, Clara Bayley


----- CLEO -----
File: divfree2d.py
Project: divfreemotion
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
Script generates input files, then runs CLEO executable "divfree2d" to create the
data to then plot for divergence free motion of superdroplets in a 2-D divergence
"""

import os
import sys
from pathlib import Path
import divfree2d_inputfiles

path2CLEO = sys.argv[1]
path2build = sys.argv[2]
configfile = sys.argv[3]

sys.path.append(path2CLEO)  # for imports from pySD package
sys.path.append(
    path2CLEO + "/examples/exampleplotting/"
)  # for imports from example plotting package

from plotssrc import pltsds, pltmoms
from pySD.sdmout_src import pyzarr, pysetuptxt, pygbxsdat

### ---------------------------------------------------------------- ###
### ----------------------- INPUT PARAMETERS ----------------------- ###
### ---------------------------------------------------------------- ###
### --- essential paths and filenames --- ###
# path and filenames for creating initial SD conditions
binpath = path2build + "/bin/"
sharepath = path2build + "/share/"
gridfile = sharepath + "df2d_dimlessGBxboundaries.dat"
initSDsfile = sharepath + "df2d_dimlessSDsinit.dat"
thermofile = sharepath + "/df2d_dimlessthermo.dat"
savefigpath = path2build + "/bin/"  # directory for saving figures

# path and file names for plotting results
setupfile = binpath + "df2d_setup.txt"
dataset = binpath + "df2d_sol.zarr"

### ---------------------------------------------------------------- ###
### ------------------- BINARY FILES GENERATION--------------------- ###
### ---------------------------------------------------------------- ###
### --- ensure build, share and bin directories exist --- ###
if path2CLEO == path2build:
    raise ValueError("build directory cannot be CLEO")
else:
    Path(path2build).mkdir(exist_ok=True)
    Path(sharepath).mkdir(exist_ok=True)
    Path(binpath).mkdir(exist_ok=True)
    Path(savefigpath).mkdir(exist_ok=True)

### --- delete any existing initial conditions --- ###
os.system("rm " + gridfile)
os.system("rm " + initSDsfile)
os.system("rm " + thermofile[:-4] + "*")

divfree2d_inputfiles.main(
    path2CLEO, path2build, configfile, gridfile, initSDsfile, thermofile
)
### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###

### ---------------------------------------------------------------- ###
### ---------------------- RUN CLEO EXECUTABLE --------------------- ###
### ---------------------------------------------------------------- ###
os.chdir(path2build)
os.system("pwd")
os.system("rm -rf " + dataset)  # delete any existing dataset
executable = path2build + "/examples/divfreemotion/src/divfree2d"
print("Executable: " + executable)
print("Config file: " + configfile)
os.system(executable + " " + configfile)
### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###

### ---------------------------------------------------------------- ###
### ------------------------- PLOT RESULTS ------------------------- ###
### ---------------------------------------------------------------- ###
# read in constants and intial setup from setup .txt file
config = pysetuptxt.get_config(setupfile, nattrs=3, isprint=True)
consts = pysetuptxt.get_consts(setupfile, isprint=True)
gbxs = pygbxsdat.get_gridboxes(gridfile, consts["COORD0"], isprint=True)

time = pyzarr.get_time(dataset)
sddata = pyzarr.get_supers(dataset, consts)
maxnsupers = pyzarr.get_totnsupers(dataset)

# 4. plot results
savename = savefigpath + "df2d_maxnsupers_validation.png"
pltmoms.plot_totnsupers(time, maxnsupers, savename=savename)

nsample = 500
savename = savefigpath + "df2d_motion2d_validation.png"
pltsds.plot_randomsample_superdrops_2dmotion(
    sddata, config["maxnsupers"], nsample, savename=savename, arrows=False
)
### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###
