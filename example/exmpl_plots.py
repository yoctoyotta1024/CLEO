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
from matplotlib.colors import LogNorm

from exmpl_figsrc import *
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

### --------------------- plot .gif animations --------------------- ###
nframes = len(time)

mom2ani = massmoms.nsupers
xlims = [0, np.amax(mom2ani)]
savename="nsupers_1dprofile"
xlabel = "number of superdroplets per gridbox"
animate1dprofile(gbxs, mom2ani, time, nframes, xlabel=xlabel, xlims=xlims,
                 color="blue", saveani=True, savedir=savedir,
                 savename=savename, fps=5)

mom2ani = massmoms.mom1 / gbxs.gbxvols[None,:,:,:]
xlims = [0, np.amax(np.mean(mom2ani, axis=(0,1)))]
savename="massconc_1dprofile"
xlabel = "Mass Concentration /g m$^{-3}$"
animate1dprofile(gbxs, mom2ani, time, nframes, xlabel=xlabel, xlims=xlims,
                 color="green", saveani=True, savedir=savedir,
                 savename=savename, fps=5)

mom2ani = np.sum(massmoms.nsupers, axis=1) # sum over y dimension
cmap="plasma"
cmapnorm = LogNorm(vmin=1, vmax=np.amax(mom2ani))
cbarlabel="number of superdroplets per gridbox"
savename="nsupers_2dcmap"
animate2dcmap(gbxs, mom2ani, time, nframes, 
                  cbarlabel=cbarlabel, cmapnorm=cmapnorm, cmap=cmap,
                  saveani=True, savedir=savedir, savename=savename, fps=5)

mom2ani = np.mean(massmoms.mom1 / gbxs.gbxvols[None,:,:,:], axis=1) # avg over y dimension
cmap="afmhot_r"
cmapnorm = LogNorm(vmin=1e-3, vmax=1e5)
savename="massconc_2dcmap"
xlabel = "Mass Concentration /g m$^{-3}$"
animate2dcmap(gbxs, mom2ani, time, nframes, 
                  cbarlabel=cbarlabel, cmapnorm=cmapnorm, cmap=cmap,
                  saveani=True, savedir=savedir, savename=savename, fps=5)                
### ---------------------------------------------------------------- ### 