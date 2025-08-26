"""
Copyright (c) 2024 MPI-M, Clara Bayley


----- CLEO -----
File: exampleplotting.py
Project: exampleplotting
Created Date: Sunday 26th November 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
"""

# %%
import os
import sys
import awkward as ak
import matplotlib.pyplot as plt
from pathlib import Path

path2cleopy = os.path.dirname(os.path.realpath(__file__)) + "/../../"
sys.path.append(path2cleopy)
from cleopy.sdmout_src import pyzarr, pysetuptxt, pygbxsdat
from plotssrc import pltsds, pltdist

# %%
### ---------------------- input parameters ------------------------ ###
### paths to data for plotting
dataset = Path("/home/m/m300950/CLEO/build/bin/SDMdata.zarr")
setuptxt = Path("/home/m/m300950/CLEO/build/bin/setup.txt")
grid_filename = Path("/home/m/m300950/CLEO/build/share/dimlessGBxboundaries.dat")

### whether and where to save figures
savefig = False
savefigpath = Path.cwd()

### individual superdroplet plotting parameters
nsample = 50

### droplet distributions plotting parameters
t2plts = [0, 100, 200, 600, 1200, 1800, 2400, 3600]
rspan = ["min", "max"]
nbins = 100

# domain mass distrib
smoothsig_mass = 0.62
perlogR_mass = True
ylog_mass = False

# domain number concentraion distrib
smoothsig_num = False
perlogR_num = False
ylog_num = True
### ---------------------------------------------------------------- ###

# %%
### ------------------------- load data ---------------------------- ###
config = pysetuptxt.get_config(setuptxt, nattrs=3, isprint=True)
consts = pysetuptxt.get_consts(setuptxt, isprint=True)
gbxs = pygbxsdat.get_gridboxes(grid_filename, consts["COORD0"], isprint=True)
ds = pyzarr.get_rawdataset(dataset)

time = pyzarr.get_time(ds)
superdrops = pyzarr.get_supers(ds, consts)
savename = ""
### ---------------------------------------------------------------- ###

# %%
### ----------------- plot individual superdroplets ---------------- ###
if savefig:
    savename = savefigpath / "randomsample_attrs.png"
pltsds.plot_randomsample_superdrops(time, superdrops, nsample, savename=savename)
plt.show()

if savefig:
    savename = savefigpath / "randomsample_2dmotion.png"
pltsds.plot_randomsample_superdrops_2dmotion(
    superdrops, nsample, savename=savename, arrows=False
)
plt.show()
superdrops.detach_time()
### ---------------------------------------------------------------- ###

# %%
### ------------------ plot droplet distributions ------------------ ###
if rspan == ["min", "max"]:
    non_nanradius = ak.nan_to_none(superdrops["radius"])
    rspan = [ak.min(non_nanradius), ak.max(non_nanradius)]

if smoothsig_mass:
    smoothsig_mass = smoothsig_mass * (config["maxnsupers"] ** (-1 / 5))
if savefig:
    savename = savefigpath / "domain_mass_distrib.png"
fig, ax = pltdist.plot_domainmass_distribs(
    time,
    superdrops,
    t2plts,
    gbxs["domainvol"],
    rspan,
    nbins,
    smoothsig=smoothsig_mass,
    perlogR=perlogR_mass,
    ylog=ylog_mass,
    savename=savename,
)
plt.show()

if smoothsig_num:
    smoothsig_num = smoothsig_num * (config["maxnsupers"] ** (-1 / 5))
if savefig:
    savename = savefigpath / "domain_numconc_distrib.png"
fig, ax = pltdist.plot_domainnumconc_distribs(
    time,
    superdrops,
    t2plts,
    gbxs["domainvol"],
    rspan,
    nbins,
    smoothsig=smoothsig_num,
    perlogR=perlogR_num,
    ylog=ylog_num,
    savename=savename,
)
plt.show()
### ---------------------------------------------------------------- ###
