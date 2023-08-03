### Author: Clara Bayley
### File: exmpl_plots.py
### python script to plot some of the output
### data from the example run of CLEO
### e.g via command line:
### python exmpl_plots.py ./build/bin/SDMdata.zarr/ ./build/bin/setup.txt ./build/share/dimlessGBxboundaries.dat 

import sys
import numpy as np
from matplotlib.colors import LogNorm, Normalize

from exmpl_figsrc import *
from exmpl_anisrc import *

### ----------------------- output data to plot -------------------- ###
dataset = sys.argv[1]
setuptxt = sys.argv[2]
gridfile = sys.argv[3]
savedir = "./build/bin/"

time = get_timemins(dataset) 
setup = Setup(setuptxt, ntime=len(time), isprint=False)
gbxs = GridBoxes(gridfile, setup, isprint=False)
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
### -------------------------------------------------------------- ###

### -------------------- plot 1-D .gif animations ------------------ ###
nframes = len(time)

def horizontal_average(data4d):
  '''avg 4-D data with dims [time, y, x, z]
  over x and y dimensions '''
  return np.mean(data4d, axis=(1,2))

mom2ani = horizontal_average(massmoms.nsupers) 
xlims = [0, np.amax(mom2ani)]
xlabel = "mean number of superdroplets per gridbox"
savename="nsupers1d"
animate1dprofile(gbxs, mom2ani, time, nframes, xlabel=xlabel, xlims=xlims,
                 color="blue", saveani=True, savedir=savedir,
                 savename=savename, fps=5)

norm = np.sum(gbxs.gbxvols, axis=0)[None,None,:,:] * 1e6 # volume [cm^3]
mom2ani = horizontal_average(massmoms.mom0/norm) 
xlims = [0, np.amax(mom2ani)]
xlabel = "mean number concentration /cm$^{-3}$"
savename="numconc1d"
animate1dprofile(gbxs, mom2ani, time, nframes, xlabel=xlabel, xlims=xlims,
                 color="green", saveani=True, savedir=savedir,
                 savename=savename, fps=5)

norm = np.sum(gbxs.gbxvols, axis=0)[None,None,:,:]  # volume [m^3]
mom2ani = horizontal_average(massmoms.mom1/norm) 
xlims = [0, np.amax(mom2ani)]
xlabel = "mean mass concentration /g m$^{-3}$"
savename="massconc1d"
animate1dprofile(gbxs, mom2ani, time, nframes, xlabel=xlabel, xlims=xlims,
                 color="green", saveani=True, savedir=savedir,
                 savename=savename, fps=5)

### -------------------------------------------------------------- ###

### -------------------- plot 2-D .gif animations ------------------ ###
nframes = len(time)

mom2ani = np.sum(massmoms.nsupers, axis=1) # sum over y dimension
cmap="plasma_r"
cmapnorm = Normalize(vmin=1, vmax=20)
cbarlabel="number of superdroplets per gridbox"
savename="nsupers2d"
animate2dcmap(gbxs, mom2ani, time, nframes, 
                  cbarlabel=cbarlabel, cmapnorm=cmapnorm, cmap=cmap,
                  saveani=True, savedir=savedir, savename=savename, fps=5)

mom2ani = np.sum(massmoms.mom0, axis=1) # sum over y dimension
norm = np.sum(gbxs.gbxvols, axis=0)[None,:,:] * 1e6 # sum over y dimension and add time dimension for broadcasting [cm^3]
mom2ani = mom2ani / norm
cmap="afmhot_r"
cmapnorm = Normalize(vmin=0, vmax=10)
cbarlabel = "number concentration /cm$^{-3}$"
savename="numconc2d"
animate2dcmap(gbxs, mom2ani, time, nframes, 
                  cbarlabel=cbarlabel, cmapnorm=cmapnorm, cmap=cmap,
                  saveani=True, savedir=savedir, savename=savename, fps=5)   

mom2ani = np.sum(massmoms.mom1, axis=1) # sum over y dimension
norm = np.sum(gbxs.gbxvols, axis=0)[None,:,:] # sum over y dimension and add time dimension for broadcasting [m^3]
mom2ani = mom2ani / norm
cmap="bone_r"
cmapnorm = LogNorm(vmin=1e-6, vmax=1e2)
cbarlabel = "mass concentration /g m$^{-3}$"
savename="massconc2d"
animate2dcmap(gbxs, mom2ani, time, nframes, 
              cbarlabel=cbarlabel, cmapnorm=cmapnorm, cmap=cmap,
              saveani=True, savedir=savedir, savename=savename, fps=5)   
             
### ---------------------------------------------------------------- ### 