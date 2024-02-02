'''
----- CLEO -----
File: rainshaft1d.py
Project: rainshaft1d
Created Date: Friday 17th November 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
Last Modified: Wednesday 17th January 2024
Modified By: CB
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
Copyright (c) 2023 MPI-M, Clara Bayley
-----
File Description:
Script compiles and runs CLEO rain1D to create the
data and plots precipitation example given constant
1-D rainshaft thermodynamics read from a file
'''
#%%
import sys
import numpy as np
import random
from pathlib import Path

path2CLEO = Path("/home/m/m301096/CLEO")
path2build = path2CLEO / "build"
configfile = path2CLEO / "examples/rainshaft1d/src/config/rain1d_config.txt"

sys.path.append(path2CLEO)  # for imports from pySD package
sys.path.append(path2CLEO / "examples/exampleplotting/") # for imports from example plotting package

from pySD.plotssrc import pltsds, pltmoms, animations
from pySD.sdmout_src import *

from pySD.gbxboundariesbinary_src import read_gbxboundaries as rgrid
from pySD.initsuperdropsbinary_src import *
from pySD.initsuperdropsbinary_src import read_initsuperdrops as rsupers
from pySD.thermobinary_src import read_thermodynamics as rthermo


from sdm_eurec4a.visulization import set_custom_rcParams

set_custom_rcParams()
# from pySD.initsuperdropsbinary_src import *

### ---------------------------------------------------------------- ###
### ----------------------- INPUT PARAMETERS ----------------------- ###
### ---------------------------------------------------------------- ###
### --- essential paths and filenames --- ###
# path and filenames for creating initial SD conditions
constsfile    = path2CLEO / "libs/cleoconstants.hpp"
binpath       = path2build / "bin/"
sharepath     = path2build / "share/"
gridfile      = sharepath / "rain1d_dimlessGBxboundaries.dat"
initSDsfile   = sharepath / "rain1d_dimlessSDsinit.dat"
thermofile    =  sharepath / "rain1d_dimlessthermo.dat"

# path and file names for plotting results
setupfile     = binpath / "rain1d_setup.txt"
dataset       = binpath / "rain1d_sol.zarr"

### --- plotting initialisation figures --- ###
isfigures   = [True, True] # booleans for [making, saving] initialisation figures
savefigpath = path2CLEO / "results/examplesolutions/no_aerosol/" # directory for saving figures
savefigpath.mkdir(exist_ok=True)

SDgbxs2plt  = list(range(39, 45))
SDgbxs2plt  = [random.choice(SDgbxs2plt)] # choose random gbx from list to plot
# %%

### ----- show (and save) plots of binary file data ----- ###
rgrid.plot_gridboxboundaries(
constsfile = str(constsfile),
gridfile   = str(gridfile),
binpath = str(savefigpath) + "/",
savefig= True
)
rthermo.plot_thermodynamics(
    constsfile = str(constsfile),
    configfile= str(configfile),
    gridfile   = str(gridfile),
    thermofile= str(thermofile),
    binpath= str(savefigpath) + "/",
    savefig= True
    )
rsupers.plot_initGBxs_distribs(
    configfile= str(configfile),
    constsfile = str(constsfile),
    initsupersfile= str(initSDsfile),
    gridfile   = str(gridfile),
    binpath= str(savefigpath) + "/",
    savefig= True,
    gbxs2plt=SDgbxs2plt
    )

# %%
### ------------------------------------------------------------ ###
### ----------------------- PLOT RESULTS ----------------------- ###
### ------------------------------------------------------------ ###
# read in constants and intial setup from setup .txt file
config = pysetuptxt.get_config(setupfile, nattrs=3, isprint=True)
consts = pysetuptxt.get_consts(setupfile, isprint=True)
gbxs = pygbxsdat.get_gridboxes(str(gridfile), consts["COORD0"], isprint=True)

time = pyzarr.get_time(str(dataset))
sddata = pyzarr.get_supers(str(dataset), consts)
totnsupers = pyzarr.get_totnsupers(str(dataset))
massmoms = pyzarr.get_massmoms(str(dataset), config["ntime"], gbxs["ndims"])

# plot figures
savename = str(savefigpath / "rain1d_totnsupers.png")
pltmoms.plot_totnsupers(time, totnsupers, savename=savename)

savename = str(savefigpath / "rain1d_domainmassmoms.png")
pltmoms.plot_domainmassmoments(time, massmoms, savename=savename)

nsample = 10
savename = str(savefigpath / "rain1d_randomsample.png")
pltsds.plot_randomsample_superdrops(time, sddata,
                                        config["totnsupers"],
                                        nsample,
                                        savename=savename)

# %%
### ----- plot 1-D .gif animations ----- ###
nframes = len(time.mins)
mom2ani = np.sum(massmoms.nsupers, axis=(1,2))
xlims = [0, np.amax(mom2ani)]
xlabel = "number of super-droplets"
savename=str(savefigpath / "rain1d_nsupers1d")
animations.animate1dprofile(gbxs, mom2ani, time.mins, nframes,
                            xlabel=xlabel, xlims=xlims,
                            color="green", saveani=True,
                            savename=savename, fps=5)

nframes = len(time.mins)
norm = gbxs["gbxvols"] * 1e6 # volume [cm^3]
mom2ani = np.sum(massmoms.mom0 / norm[None,:], axis=(1,2))
xlims = [0, np.amax(mom2ani)]
xlabel = "number concentration /cm$^{-3}$"
savename=str(savefigpath / "rain1d_numconc1d")
animations.animate1dprofile(gbxs, mom2ani, time.mins, nframes,
                            xlabel=xlabel, xlims=xlims,
                            color="green", saveani=True,
                            savename=savename, fps=5)

nframes = len(time.mins)
norm = gbxs["gbxvols"] # volume [m^3]
mom2ani = np.sum(massmoms.mom1/ norm[None,:], axis=(1,2))
xlims = [0, np.amax(mom2ani)]
xlabel = "mass concentration /g m$^{-3}$"
savename=str(savefigpath / "rain1d_massconc1d")
animations.animate1dprofile(gbxs, mom2ani, time.mins, nframes,
                            xlabel=xlabel, xlims=xlims,
                            color="green", saveani=True,
                            savename=savename, fps=5)

print(savefigpath)
### ------------------------------------------------------------ ###
### ------------------------------------------------------------ ###

# %%
