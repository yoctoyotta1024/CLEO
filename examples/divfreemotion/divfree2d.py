'''
----- CLEO -----
File: divfree2d.py
Project: divfreemotion
Created Date: Friday 17th November 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
Last Modified: Wednesday 22nd November 2023
Modified By: CB
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
Copyright (c) 2023 MPI-M, Clara Bayley
-----
File Description:
Script compiles and runs CLEO divfree2D to create the
data and plots for divergence free motion of
superdroplets in a 2-D divergence free wind field
'''

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

path2CLEO = sys.argv[1]
path2build = sys.argv[2]
configfile = sys.argv[3]

sys.path.append(path2CLEO)  # for imports from pySD package
sys.path.append(path2CLEO+"/examples/exampleplotting/") # for imports from example plotting package

from plotssrc import pltsds, pltmoms
from pySD.sdmout_src import *
from pySD.initsuperdropsbinary_src import *
from pySD.gbxboundariesbinary_src import read_gbxboundaries as rgrid
from pySD.gbxboundariesbinary_src import create_gbxboundaries as cgrid
from pySD.initsuperdropsbinary_src import initattributes as iattrs
from pySD.initsuperdropsbinary_src import radiiprobdistribs as rprobs 
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
gridfile = sharepath+"df2d_dimlessGBxboundaries.dat"
initSDsfile = sharepath+"df2d_dimlessSDsinit.dat"
thermofile =  sharepath+"/df2d_dimlessthermo.dat"

# path and file names for plotting results
setupfile = binpath+"df2d_setup.txt"
dataset = binpath+"df2d_sol.zarr"

### --- plotting initialisation figures --- ###
isfigures = [True, True] # booleans for [making, saving] initialisation figures
savefigpath = path2build+"/bin/" # directory for saving figures
SDgbxs2plt = [0] # gbxindex of SDs to plot (nb. "all" can be very slow)

### --- settings for 2-D gridbox boundaries --- ###
zgrid = [0, 1500, 50]     # evenly spaced zhalf coords [zmin, zmax, zdelta] [m]
xgrid = [0, 1500, 50]     # evenly spaced xhalf coords [m]
ygrid = np.array([0, 20])  # array of yhalf coords [m]

### --- settings for initial superdroplets --- ###
# settings for initial superdroplet coordinates
zlim = 750        # max z coord of superdroplets
npergbx = 8       # number of superdroplets per gridbox 

# [min, max] range of initial superdroplet radii (and implicitly solute masses)
rspan                = [3e-9, 3e-6] # [m]

# settings for initial superdroplet multiplicies
# (from bimodal Lognormal distribution)
geomeans             = [0.02e-6, 0.15e-6]               
geosigs              = [1.4, 1.6]                    
scalefacs            = [6e6, 4e6]   
numconc = np.sum(scalefacs)

### --- settings for 2D Thermodyanmics --- ###
PRESS0 = 100000 # [Pa]
THETA = 298.15  # [K]
qcond = 0.0     # [Kg/Kg]
WMAX = 0.6      # [m/s]
VVEL = None     # [m/s]
Zlength = 1500  # [m]
Xlength = 1500  # [m]
qvapmethod = "sratio"
Zbase = 750     # [m]
sratios = [1.0, 1.0] # s_ratio [below, above] Zbase
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
  if isfigures[1]:
    Path(savefigpath).mkdir(exist_ok=True) 
os.system("rm "+gridfile)
os.system("rm "+initSDsfile)
os.system("rm "+thermofile[:-4]+"*")

### ----- write gridbox boundaries binary ----- ###
cgrid.write_gridboxboundaries_binary(gridfile, zgrid, xgrid, ygrid, constsfile)
rgrid.print_domain_info(constsfile, gridfile)

### ----- write thermodyanmics binaries ----- ###
thermodyngen = thermogen.SimpleThermo2Dflowfield(configfile, constsfile, PRESS0,
                                                THETA, qvapmethod, sratios, Zbase,
                                                qcond, WMAX, Zlength, Xlength,
                                                VVEL)
cthermo.write_thermodynamics_binary(thermofile, thermodyngen, configfile,
                                    constsfile, gridfile)


### ----- write initial superdroplets binary ----- ###
nsupers = iattrs.nsupers_at_domain_base(gridfile, constsfile, npergbx, zlim)
coord3gen = iattrs.SampleCoordGen(True) # sample coord3 randomly
coord1gen = iattrs.SampleCoordGen(True) # sample coord1 randomly
coord2gen = None                        # do not generate superdroplet coord2s
radiiprobdist = rprobs.LnNormal(geomeans, geosigs, scalefacs)
radiigen = iattrs.SampleDryradiiGen(rspan, True) # randomly sample radii from rspan [m]

initattrsgen = iattrs.InitManyAttrsGen(radiigen, radiiprobdist,
                                        coord3gen, coord1gen, coord2gen)
csupers.write_initsuperdrops_binary(initSDsfile, initattrsgen, 
                                      configfile, constsfile,
                                      gridfile, nsupers, numconc)

### ----- show (and save) plots of binary file data ----- ###
if isfigures[0]:
  rgrid.plot_gridboxboundaries(constsfile, gridfile,
                               savefigpath, isfigures[1])
  rthermo.plot_thermodynamics(constsfile, configfile, gridfile,
                              thermofile, savefigpath, isfigures[1])
  rsupers.plot_initGBxsdistribs(configfile, constsfile, initSDsfile,
                              gridfile, savefigpath, isfigures[1],
                              SDgbxs2plt) 
  plt.close()
### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###

### ---------------------------------------------------------------- ###
### -------------------- COMPILE AND RUN CLEO ---------------------- ###
### ---------------------------------------------------------------- ###
os.chdir(path2build)
os.system('pwd')
os.system('rm -rf '+dataset)
os.system('make clean && make -j 64 divfree2D')
executable = path2build+'/examples/divfreemotion/src/divfree2D'
os.system(executable + ' ' + configfile)
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
totnsupers =pyzarr.get_totnsupers(dataset)

# 4. plot results
savename = savefigpath + "df2d_nsupers_validation.png"
pltmoms.plot_totnsupers(time, totnsupers, savename=savename)

nsample = 500
savename = savefigpath + "df2d_motion2d_validation.png"
pltsds.plot_randomsample_superdrops_2dmotion(sddata,
                                                 config["totnsupers"],
                                                 nsample,
                                                 savename=savename,
                                                 arrows=False)
### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###
