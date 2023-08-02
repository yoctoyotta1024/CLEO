### Author: Clara Bayley
### File: exmpl_plots.py
### python script to plot some of the output
### data from the example run of CLEO
### e.g via command line:
### python exmpl_plots.py ./build/bin/SDMdata.zarr/ ./build/bin/setup.txt ./build/share/dimlessGBxboundaries.dat 

import sys
import numpy as np
import xarray as xr
import awkward as ak
import matplotlib.pyplot as plt

from exmpl_plotssrc import *
from exmpl_anisrc import *

### ----------------------- output data to plot -------------------- ###
dataset = sys.argv[1]
setuptxt = sys.argv[2]
gridfile = sys.argv[3]
savedir = "./build/bin/"

setup = Setup(setuptxt, isprint=False)
gbxs = GridBoxes(gridfile, setup, isprint=False)
time = get_timemins(dataset) 
massmoms = MassMoments(dataset, setup, gbxs.ndims)
sddata = SDData(dataset)
precip = SurfPrecip(dataset, setup, gbxs)
### ---------------------------------------------------------------- ###

### ----------------------- plot .png figures ---------------------- ###
savename = "domainmassmoments.png"
argsdict = {'time': time,
            'massmoms': massmoms}
genericplot(plot_domainmassmoments, argsdict,
            figsize=(10,6), savefig=True,
            savedir=savedir, savename=savename)

savename = "randomsampleSDs.png"
argsdict = {'time': time,
            'sddata': sddata,
            'nsample':  50}
genericplot(plot_randomsampleSDs, argsdict,
            figsize=(10,6), savefig=True,
            savedir=savedir, savename=savename)

savename = "surfaceprecip.png"
argsdict = {'time': time,
            'precip': precip}
genericplot(plot_surfaceprecip, argsdict,
            figsize=(10,6), savefig=True,
            savedir=savedir, savename=savename)
### ---------------------------------------------------------------- ###

### --------------------- plot .mp4 animations --------------------- ###



### ---------------------------------------------------------------- ###