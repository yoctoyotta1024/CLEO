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
import yaml

import xarray as xr
import matplotlib.pyplot as plt

path2CLEO = Path("/home/m/m301096/CLEO")
path2build = path2CLEO / "build"
configfile = path2CLEO / "examples/eurec4a_rainshaft1d/src/config/rain1d_config.txt"

yaml_config_file = Path("/home/m/m301096/repositories/sdm-eurec4a/data/model/input/example_input.yaml")
with open(yaml_config_file, 'r') as f:
    config_yaml = yaml.safe_load(f)



sys.path.append(path2CLEO)  # for imports from pySD package
sys.path.append(path2CLEO / "examples/exampleplotting") # for imports from example plotting package

# from plotssrc import pltsds, pltmoms, animations
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
cloud_id = config_yaml['cloud']['cloud_id']
savefigpath = path2CLEO / "results/plots/yaml_config/" / f"cloud_{cloud_id}" # directory for saving figures
savefigpath.mkdir(exist_ok=True, parents=True)

SDgbxs2plt  = list(range(39, 45))
SDgbxs2plt  = [random.choice(SDgbxs2plt)] # choose random gbx from list to plot
# %%

# read in constants and intial setup from setup .txt file
config = pysetuptxt.get_config(setupfile, nattrs=3, isprint=True)
consts = pysetuptxt.get_consts(setupfile, isprint=True)
gbxs = pygbxsdat.get_gridboxes(str(gridfile), consts["COORD0"], isprint=True)

time = pyzarr.get_time(str(dataset))
sddata = pyzarr.get_supers(str(dataset), consts)
totnsupers = pyzarr.get_totnsupers(str(dataset))
massmoms = pyzarr.get_massmoms(str(dataset), config["ntime"], gbxs["ndims"])



# simple_ds = xr.open_dataset(dataset, engine="zarr",
#                                 consolidated=False)

# # %%
# # To extract the data of all superdroplets, one needs to use the from pySD.sdmout_src import sdtracing

# from pySD.sdmout_src import sdtracing
# all_sd_ids = np.arange(0, 10) # int(config["totnsupers"]))

# attributes = ["xi", "radius", "coord3"]

# list_dataarrays = []
# for attr in attributes:
#     data = sdtracing.attribute_for_superdroplets_sample(
#         sddata,
#         attr,
#         ids=list(all_sd_ids)
#     )
#     da = xr.DataArray(
#         data=data,
#         dims=["time", "sd_id"],
#         coords={"time": simple_ds.time, "sd_id": all_sd_ids}
#     )
#     da.name = attr
#     list_dataarrays.append(da)


# # %%
# ds = xr.merge(list_dataarrays)

# %%
### ------------------------------------------------------------ ###
### ----------------------- PLOT RESULTS ----------------------- ###
### ------------------------------------------------------------ ###

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
