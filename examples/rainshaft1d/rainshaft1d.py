'''
----- CLEO -----
File: rainshaft1d.py
Project: rainshaft1d
Created Date: Friday 17th November 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
Last Modified: Tuesday 9th January 2024
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
gridfile = sharepath+"rain1d_dimlessGBxboundaries.dat"
initSDsfile = sharepath+"rain1d_dimlessSDsinit.dat"
thermofile =  sharepath+"rain12d_dimlessthermo.dat"

# path and file names for plotting results
setupfile = binpath+"rain1d_setup.txt"
dataset = binpath+"rain1d_sol.zarr"

### --- plotting initialisation figures --- ###
isfigures = [True, True] # booleans for [making, saving] initialisation figures
savefigpath = path2build+"/bin/" # directory for saving figures
SDgbxs2plt = [0] # gbxindex of SDs to plot (nb. "all" can be very slow)

### --- settings for 1-D gridbox boundaries --- ###
zgrid = [0, 1500, 20]      # evenly spaced zhalf coords [zmin, zmax, zdelta] [m]
xgrid = np.array([0, 20])  # array of xhalf coords [m]
ygrid = np.array([0, 20])  # array of yhalf coords [m]

### --- settings for 1-D Thermodynamics --- ###
PRESS0 = 101315 # [Pa]
THETA = 298.15  # [K]
qcond = 0.0     # [Kg/Kg]
WMAX = 0.6      # [m/s]
VVEL = None     # [m/s]
Zlength = 1500  # [m]
Xlength = 1500  # [m]
qvapmethod = "sratio"
Zbase = 750 # [m]
sratios = [0.99, 1.0025] # s_ratio [below, above] Zbase
moistlayer=False

### --- settings for initial superdroplets --- ###
# settings for initial superdroplet coordinates
zlim = 1000       # min z coord of superdroplets
npergbx = 256    # number of superdroplets per gridbox 

# [min, max] range of initial superdroplet radii (and implicitly solute masses)
rspan                = [1e-8, 1e-4]                  # random sample of radii in this range [m]
monodryr             = 1e-8                          # all SDs have this same dryradius [m]

# settings for initial superdroplet multiplicies
reff                 = 7e-6                     # effective radius [m]
nueff                = 0.08                     # effective variance 
rdist1 = probdists.ClouddropsHansenGamma(reff, nueff)

nrain                = 3000                         # raindrop concentration [m^-3]
qrain                = 0.9                          # rainwater content [g/m^3]
dvol                 = 8e-4                         # mean volume diameter [m]
rdist2 = probdists.RaindropsGeoffroyGamma(nrain, qrain, dvol)

distribs = [rdist1, rdist2]
scalefacs = [1000, 1]                               # relative abundance of 2 distributions
numconc = 1e9                                       # [m^3]
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
os.system("rm "+gridfile)
os.system("rm "+initSDsfile)
os.system("rm "+thermofile[:-4]+"*")

### ----- write gridbox boundaries binary ----- ###
cgrid.write_gridboxboundaries_binary(gridfile, zgrid, xgrid, ygrid, constsfile)
rgrid.print_domain_info(constsfile, gridfile)

### ----- write thermodynamics binaries ----- ###
thermodyngen = thermogen.ConstHydrostaticAdiabat(configfile, constsfile, PRESS0, 
                                                 THETA, qvapmethod, sratios, Zbase,
                                                 qcond, WMAX, Zlength, Xlength,
                                                 VVEL, moistlayer)
cthermo.write_thermodynamics_binary(thermofile, thermodyngen, configfile,
                                    constsfile, gridfile)


### ----- write initial superdroplets binary ----- ###
nsupers = crdgens.nsupers_at_domain_top(gridfile, constsfile, npergbx, zlim)
coord3gen = crdgens.SampleCoordGen(True) # sample coord3 randomly
coord1gen = None                        # do not generate superdroplet coord2s
coord2gen = None                        # do not generate superdroplet coord2s

xiprobdist = probdists.CombinedRadiiProbDistribs(distribs, scalefacs)
radiigen = rgens.SampleLog10RadiiGen(rspan) # randomly sample radii from rspan [m]
dryradiigen  =  rgens.MonoAttrGen(monodryr)             # all SDs have the same dryradius [m]

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
### -------------------- COMPILE AND RUN CLEO ---------------------- ###
### ---------------------------------------------------------------- ###
# 2. compile and the run model
os.chdir(path2build)
os.system('pwd')
os.system('rm -rf '+dataset)
os.system('make clean && make -j 64 rain1D')
executable = path2build+'/examples/constthermo2d/src/rain1D'
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
savename = savefigpath + "rain1d_totnsupers.png"
pltmoms.plot_totnsupers(time, totnsupers, savename=savename)

savename = savefigpath + "rain1d_domainmassmoms.png"
pltmoms.plot_domainmassmoments(time, massmoms, savename=savename)

nsample = 500
savename = savefigpath + "rain1d_randomsample.png"
pltsds.plot_randomsample_superdrops(time, sddata,
                                        config["totnsupers"],
                                        nsample,
                                        savename=savename)

### ----- plot 1-D .gif animations ----- ###
nframes = len(time.mins)
mom2ani = np.sum(massmoms.nsupers, axis=(1,2))
xlims = [0, np.amax(mom2ani)]
xlabel = "number of super-droplets"
savename=savefigpath+"rain1d_nsupers1d"
animations.animate1dprofile(gbxs, mom2ani, time.mins, nframes,
                            xlabel=xlabel, xlims=xlims,
                            color="green", saveani=True,
                            savename=savename, fps=5)   

nframes = len(time.mins)
norm = gbxs["gbxvols"] * 1e6 # volume [cm^3]
mom2ani = np.sum(massmoms.mom0, axis=(1,2)) / norm[None,:]
xlims = [0, np.amax(mom2ani)]
xlabel = "number concentration /cm$^{-3}$"
savename=savefigpath+"rain1d_numconc1d"
animations.animate1dprofile(gbxs, mom2ani, time.mins, nframes,
                            xlabel=xlabel, xlims=xlims,
                            color="green", saveani=True,
                            savename=savename, fps=5)

nframes = len(time.mins)
norm = gbxs["gbxvols"] # volume [m^3]
mom2ani = np.sum(massmoms.mom1, axis=(1,2)) / norm[None,:]
xlims = [0, np.amax(mom2ani)]
xlabel = "mass concentration /g m$^{-3}$"
savename=savefigpath+"rain1d_massconc1d"
animations.animate1dprofile(gbxs, mom2ani, time.mins, nframes,
                            xlabel=xlabel, xlims=xlims,
                            color="green", saveani=True,
                            savename=savename, fps=5)                        
### ------------------------------------------------------------ ###
### ------------------------------------------------------------ ###                                