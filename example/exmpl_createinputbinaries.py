### Author: Clara Bayley
### File: exmpl_createinputbinaries.py
### python script to generate some input
### binary files for an example of runing CLEO

import sys
import numpy as np
from pathlib import Path
sys.path.append("..") # used to include directory containing pySD package in python interpreter PATH

from pySD.gbxboundariesbinary_src import create_gbxboundaries as  cgrid
from pySD.gbxboundariesbinary_src import read_gbxboundaries as rgrid
from pySD.thermobinary_src import thermogen
from pySD.thermobinary_src import create_thermodynamics as cthermo
from pySD.thermobinary_src import read_thermodynamics as rthermo
from pySD.initsuperdropsbinary_src import initattributes as iSDs
from pySD.initsuperdropsbinary_src import radiiprobdistribs as rprobs 
from pySD.initsuperdropsbinary_src import create_initsuperdrops as csupers 
from pySD.initsuperdropsbinary_src import read_initsuperdrops as rsupers 

### ---------------------------------------------------------------- ###
### ----------------------- INPUT PARAMETERS ----------------------- ###
### ---------------------------------------------------------------- ###

path2CLEO = sys.argv[1]  # absolute or relative path of CLEO directory
path2build = sys.argv[2] # absolute or relative path to build directory
configfile = sys.argv[3] # absolute or relative path to configuration file
isfigures = [True, True] # booleans for [making+showing, saving] figures
SDgbxs2plt = [0] # indexes of GBx index of SDs to plot (nb. "all" can be very slow)

### --- essential paths and filenames --- ###
# where to find constants and config files
constsfile = path2CLEO+"/libs/claras_SDconstants.hpp"

# where to write input files (note these should match configfile)
path2share = path2build+"/share"
gridfile = path2share+"/dimlessGBxboundaries.dat"
thermofile =  path2share+"/dimlessthermo.dat"
initSDsfile = path2share+"/dimlessSDsinit.dat"

# directory to save figures of initial conditions in
savefigpath = path2build+"/bin/"


### --- settings for gridbox boundaries --- ###
zgrid = [0, 1500, 100]      # evenly spaced zhalf coords [zmin, zmax, zdelta] [m]
xgrid = [0, 1500, 100]     # evenly spaced xhalf coords [m]
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
PRESS0 = 101500 # [Pa]
THETA = 289 # [K]
qcond = 0.0 # [Kg/Kg]
WMAX = 0.6 # [m/s]
VVEL = None # [m/s]
Zlength = 1500 # [m]
Xlength = 1500 # [m]
qvapmethod = "sratio"
Zbase = 750 # [m]
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
  Path(path2build+"/share/").mkdir(exist_ok=True) 
  Path(path2build+"/bin/").mkdir(exist_ok=True) 


### ----- write gridbox boundaries binary ----- ###
cgrid.write_gridboxboundaries_binary(gridfile, zgrid, xgrid, ygrid, constsfile)
rgrid.print_domain_info(constsfile, gridfile)


### ----- write thermodyanmics binaries ----- ###
thermodyngen = thermogen.ConstHydrostaticAdiabat(configfile, constsfile, PRESS0, 
                                        THETA, qvapmethod, sratios, Zbase,
                                        qcond, WMAX, Zlength, Xlength,
                                        VVEL, moistlayer)
cthermo.write_thermodynamics_binary(thermofile, thermodyngen, configfile,
                                    constsfile, gridfile)


### ----- write initial superdroplets binary ----- ###
nsupers = iSDs.nsupers_at_domain_base(gridfile, constsfile, npergbx, zlim)
coord3gen = iSDs.SampleCoordGen(True) # sample coord3 randomly
coord1gen = iSDs.SampleCoordGen(True) # sample coord1 randomly
coord2gen = None                        # do not generate superdroplet coord2s
radiiprobdist = rprobs.LnNormal(geomeans, geosigs, scalefacs)
radiigen = iSDs.SampleDryradiiGen(rspan, True) # randomly sample radii from rspan [m]

initattrsgen = iSDs.InitManyAttrsGen(radiigen, radiiprobdist,
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
  rsupers.plot_initGBxsdistribs(configfile, constsfile, initSDsfile,
                              gridfile, savefigpath, isfigures[1],
                              SDgbxs2plt) 
### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###