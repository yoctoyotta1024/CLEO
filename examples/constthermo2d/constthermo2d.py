'''
----- CLEO -----
File: constthermo2d.py
Project: constthermo2d
Created Date: Friday 17th November 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
Last Modified: Friday 3rd May 2024
Modified By: CB
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
Copyright (c) 2023 MPI-M, Clara Bayley
-----
File Description:
Script generatees input files, runs CLEO executable "const2D" to create
data and then plots precipitation example given 2-D flow field and
constant thermodynamics read from a file.
'''

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib.colors import LogNorm, Normalize

path2CLEO = sys.argv[1]
path2build = sys.argv[2]
configfile = sys.argv[3]

sys.path.append(path2CLEO)  # for imports from pySD package
sys.path.append(path2CLEO+"/examples/exampleplotting/") # for imports from example plotting package

from plotssrc import pltsds, pltmoms, animations
from pySD.sdmout_src import *
from pySD.gbxboundariesbinary_src import read_gbxboundaries as rgrid
from pySD.gbxboundariesbinary_src import create_gbxboundaries as cgrid
from pySD.initsuperdropsbinary_src import *
from pySD.initsuperdropsbinary_src import create_initsuperdrops as csupers
from pySD.initsuperdropsbinary_src import read_initsuperdrops as rsupers
from pySD.thermobinary_src import thermogen
from pySD.thermobinary_src import create_thermodynamics as cthermo
from pySD.thermobinary_src import read_thermodynamics as rthermo

### ---------------------------------------------------------------- ###
### ----------------------- INPUT PARAMETERS ----------------------- ###
### ---------------------------------------------------------------- ###
### --- essential paths and filenames --- ###
# path and filenames for creating initial SD conditions
constsfile = path2CLEO+"/libs/cleoconstants.hpp"
binpath = path2build+"/bin/"
sharepath = path2build+"/share/"
gridfile = sharepath+"const2d_dimlessGBxboundaries.dat"
initSDsfile = sharepath+"const2d_dimlessSDsinit.dat"
thermofile =  sharepath+"/const2d_dimlessthermo.dat"

# path and file names for plotting results
setupfile = binpath+"const2d_setup.txt"
dataset = binpath+"const2d_sol.zarr"

### --- plotting initialisation figures --- ###
isfigures = [True, True] # booleans for [making, saving] initialisation figures
savefigpath = path2build+"/bin/" # directory for saving figures
SDgbxs2plt = [0] # gbxindex of SDs to plot (nb. "all" can be very slow)

### --- settings for 2-D gridbox boundaries --- ###
zgrid = [0, 1500, 75]      # evenly spaced zhalf coords [zmin, zmax, zdelta] [m]
xgrid = [0, 1500, 75]      # evenly spaced xhalf coords [m]
ygrid = np.array([0, 20])  # array of yhalf coords [m]

### --- settings for initial superdroplets --- ###
# settings for initial superdroplet coordinates
zlim = 500        # max z coord of superdroplets
npergbx = 8       # number of superdroplets per gridbox

# [min, max] range of initial superdroplet radii (and implicitly solute masses)
rspan                = [3e-9, 3e-6] # [m]

# settings for initial superdroplet multiplicies
# (from bimodal Lognormal distribution)
geomeans             = [0.02e-6, 0.15e-6]
geosigs              = [1.4, 1.6]
scalefacs            = [6e6, 4e6]
numconc = np.sum(scalefacs)

### --- settings for 2D Thermodynamics --- ###
PRESS0 = 101315     # [Pa]
THETA = 288.15      # [K]
qcond = 0.0         # [Kg/Kg]
WMAX = 0.6          # [m/s]
VVEL = None         # [m/s]
Zlength = 1500      # [m]
Xlength = 1500      # [m]
qvapmethod = "sratio"
Zbase = 750         # [m]
sratios = [0.99, 1.0025] # s_ratio [below, above] Zbase
moistlayer=False
### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###


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

### --- delete any existing initial conditions --- ###
os.system("rm "+gridfile)
os.system("rm "+initSDsfile)
os.system("rm "+thermofile[:-4]+"*")

### ----- write gridbox boundaries binary ----- ###
cgrid.write_gridboxboundaries_binary(gridfile, zgrid, xgrid, ygrid, constsfile)
rgrid.print_domain_info(constsfile, gridfile)

### ----- write thermodynamics binaries ----- ###
thermodyngen = thermogen.ConstDryHydrostaticAdiabat(configfile, constsfile, PRESS0,
                                                    THETA, qvapmethod, sratios, Zbase,
                                                    qcond, WMAX, Zlength, Xlength,
                                                    VVEL, moistlayer)
cthermo.write_thermodynamics_binary(thermofile, thermodyngen, configfile,
                                    constsfile, gridfile)


### ----- write initial superdroplets binary ----- ###
nsupers = crdgens.nsupers_at_domain_base(gridfile, constsfile, npergbx, zlim)
coord3gen = crdgens.SampleCoordGen(True) # sample coord3 randomly
coord1gen = crdgens.SampleCoordGen(True) # sample coord1 randomly
coord2gen = None                        # do not generate superdroplet coord2s
xiprobdist = probdists.LnNormal(geomeans, geosigs, scalefacs)
radiigen = rgens.SampleLog10RadiiGen(rspan) # randomly sample radii from rspan [m]
dryradiigen = dryrgens.ScaledRadiiGen(1.0)

initattrsgen = attrsgen.AttrsGenerator(radiigen, dryradiigen, xiprobdist,
                                        coord3gen, coord1gen, coord2gen)
csupers.write_initsuperdrops_binary(initSDsfile, initattrsgen,
                                      configfile, constsfile,
                                      gridfile, nsupers, numconc)

### ----- show (and save) plots of binary file data ----- ###
if isfigures[0]:
  if isfigures[1]:
    Path(savefigpath).mkdir(exist_ok=True)
  rgrid.plot_gridboxboundaries(constsfile, gridfile,
                               savefigpath, isfigures[1])
  rthermo.plot_thermodynamics(constsfile, configfile, gridfile,
                              thermofile, savefigpath, isfigures[1])
  rsupers.plot_initGBxs_distribs(configfile, constsfile, initSDsfile,
                              gridfile, savefigpath, isfigures[1],
                              SDgbxs2plt)
### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###

### ---------------------------------------------------------------- ###
### -------------------- RUN CLEO EXECUTABLE ------------------- ###
### ---------------------------------------------------------------- ###
os.chdir(path2build)
os.system('pwd')
os.system('rm -rf '+dataset) # delete any existing dataset
executable = path2build+'/examples/constthermo2d/src/const2D'
print('Executable: '+executable)
print('Config file: '+configfile)
os.system(executable + ' ' + configfile)
### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###

### ------------------------------------------------------------ ###
### ----------------------- PLOT RESULTS ----------------------- ###
### ------------------------------------------------------------ ###
# read in constants and intial setup from setup .txt file
config = pysetuptxt.get_config(setupfile, nattrs=3, isprint=True)
consts = pysetuptxt.get_consts(setupfile, isprint=True)
gbxs = pygbxsdat.get_gridboxes(gridfile, consts["COORD0"], isprint=True)

time = pyzarr.get_time(dataset)
sddata = pyzarr.get_supers(dataset, consts)
totnsupers = pyzarr.get_totnsupers(dataset)
massmoms = pyzarr.get_massmoms(dataset, config["ntime"], gbxs["ndims"])

# plot figures
savename = savefigpath + "const2d_totnsupers.png"
pltmoms.plot_totnsupers(time, totnsupers, savename=savename)

savename = savefigpath + "const2d_domainmassmoms.png"
pltmoms.plot_domainmassmoments(time, massmoms, savename=savename)

nsample = 500
savename = savefigpath + "const2d_randomsample.png"
pltsds.plot_randomsample_superdrops(time, sddata,
                                        config["maxnsupers"],
                                        nsample,
                                        savename=savename)

savename = savefigpath + "const2d_motion2d.png"
pltsds.plot_randomsample_superdrops_2dmotion(sddata,
                                                 config["maxnsupers"],
                                                 nsample,
                                                 savename=savename,
                                                 arrows=False)

### ----- plot 1-D .gif animations ----- ###
def horizontal_average(data4d):
  '''avg 4-D data with dims [time, y, x, z]
  over x and y dimensions '''
  return np.mean(data4d, axis=(1,2))

nframes = len(time.mins)
norm = np.sum(gbxs["gbxvols"], axis=0)[None,None,:,:] * 1e6 # volume [cm^3]
mom2ani = horizontal_average(massmoms.mom0/norm)
xlims = [0, np.amax(mom2ani)]
xlabel = "mean number concentration /cm$^{-3}$"
savename=savefigpath+"const2d_numconc1d"
animations.animate1dprofile(gbxs, mom2ani, time.mins, nframes,
                            xlabel=xlabel, xlims=xlims,
                            color="green", saveani=True,
                            savename=savename, fps=5)

### ----- plot 2-D .gif animations ----- ###
nframes = len(time.mins)
mom2ani = np.sum(massmoms.nsupers, axis=1) # sum over y dimension
cmap="plasma_r"
cmapnorm = Normalize(vmin=1, vmax=20)
cbarlabel="number of superdroplets per gridbox"
savename=savefigpath+"const2d_nsupers2d"
animations.animate2dcmap(gbxs, mom2ani, time.mins, nframes,
                  cbarlabel=cbarlabel, cmapnorm=cmapnorm, cmap=cmap,
                  saveani=True, savename=savename, fps=5)

nframes = len(time.mins)
mom2ani = np.sum(massmoms.mom1, axis=1) # sum over y dimension
norm = np.sum(gbxs["gbxvols"], axis=0)[None,:,:] # sum over y dimension and add time dimension for broadcasting [m^3]
mom2ani = mom2ani / norm
cmap="bone_r"
cmapnorm = LogNorm(vmin=1e-6, vmax=1e2)
cbarlabel = "mass concentration /g m$^{-3}$"
savename=savefigpath+"const2d_massconc2d"
animations.animate2dcmap(gbxs, mom2ani, time.mins, nframes,
              cbarlabel=cbarlabel, cmapnorm=cmapnorm, cmap=cmap,
              saveani=True, savename=savename, fps=5)
### ------------------------------------------------------------ ###
### ------------------------------------------------------------ ###
