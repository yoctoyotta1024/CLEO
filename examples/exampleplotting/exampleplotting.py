'''
----- CLEO -----
File: exampleplotting.py
Project: exampleplotting
Created Date: Sunday 26th November 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
Last Modified: Sunday 26th November 2023
Modified By: CB
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
Copyright (c) 2023 MPI-M, Clara Bayley
-----
File Description:
'''


import sys
import numpy as np
import xarray as xr
import awkward as ak
import matplotlib.pyplot as plt

sys.path.append('../..//CLEO/') # path to pySD
from pySD.sdmout_src import *
from pySD.sdmout_src import sdtracing
from plotssrc import pltmoms, pltsds, pltdist

### ---------------------- input parameters ------------------------ ###
### paths to data for plotting
dataset = "/home/m/m300950/CLEO/build/bin/SDMdata.zarr"
setuptxt = "/home/m/m300950/CLEO/build/bin/setup.txt"
gridfile = "/home/m/m300950/CLEO/build/share/dimlessGBxboundaries.dat"

### whether and where to save figures
savefig = True
savefigpath = "./"

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

### ------------------------- load data ---------------------------- ###
config = pysetuptxt.get_config(setuptxt, nattrs=3, isprint=True)
consts = pysetuptxt.get_consts(setuptxt, isprint=True)
gbxs = pygbxsdat.get_gridboxes(gridfile, consts["COORD0"], isprint=True)
ds = pyzarr.get_rawdataset(dataset)

time = pyzarr.get_time(ds)
thermodata = pyzarr.get_thermodata(ds, config["ntime"], gbxs["ndims"], consts)
sddata = pyzarr.get_supers(ds, consts)
gbxindex = pyzarr.get_gbxindex(ds, gbxs["ndims"])
totnsupers = pyzarr.get_totnsupers(ds)
nsupers = pyzarr.get_nsupers(ds, config["ntime"], gbxs["ndims"])
### ---------------------------------------------------------------- ###

### ----------------- plot individual superdroplets ---------------- ###
savename = ""
if savefig:
  savename = savefigpath+"/randomsample_attrs.png"
pltsds.plot_randomsample_superdrops(time, sddata, totnsupers[0], 
                                    nsample, savename=savename)
if savefig:
  savename = savefigpath+"/randomsample_2dmotion.png"
pltsds.plot_randomsample_superdrops_2dmotion(sddata, totnsupers[0], nsample, arrows=False)
### ---------------------------------------------------------------- ###

### ------------------ plot droplet distributions ------------------ ###
if rspan == ["min", "max"]:
  rspan = [np.nanmin(sddata["radius"]), np.nanmax(sddata["radius"])]

if smoothsig_mass:
  smoothsig_mass = smoothsig_mass*(config["totnsupers"]**(-1/5))
if savefig:
  savename = savefigpath+"/domain_mass_distrib.png"
fig, ax = pltdist.plot_domainmass_distribs(time.secs, sddata, t2plts, 
                                     gbxs["domainvol"], rspan, nbins,
                                     smoothsig=smoothsig_mass,
                                     perlogR=perlogR_mass,
                                     ylog=ylog_mass,
                                     savename=savename)

if smoothsig_num:
  smoothsig_num = smoothsig_num*(config["totnsupers"]**(-1/5))
if savefig:
  savename = savefigpath+"/domain_numconc_distrib.png"
fig, ax = pltdist.plot_domainnumconc_distribs(time.secs, sddata, t2plts, 
                                     gbxs["domainvol"], rspan, nbins,
                                     smoothsig=smoothsig_num,
                                     perlogR=perlogR_num,
                                     ylog=ylog_num,
                                     savename=savename)
### ---------------------------------------------------------------- ###