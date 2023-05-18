### Author: Clara Bayley
### File: quickcreate_inputbinaries.py
### python script to generate some input
### binary files for an example of runing CLEO

import sys
import numpy as np
from pathlib import Path

from pySD.gbxboundariesbinary_src import create_gbxboundaries as  cgrid
from pySD.gbxboundariesbinary_src import read_gbxboundaries as rgrid
from pySD.thermobinary_src import thermogen as gthermo
from pySD.thermobinary_src import create_thermodynamics as cthermo
from pySD.thermobinary_src import read_thermodynamics as rthermo
from pySD.initsuperdropsbinary_src import initattributes as iattrs
from pySD.initsuperdropsbinary_src import radiiprobdistribs as rprobs 
from pySD.initsuperdropsbinary_src import create_initsuperdrops as csupers 
from pySD.initsuperdropsbinary_src import read_initsuperdrops as rsupers 

### ---------------------------------------------------------------- ###
### ----------------------- INPUT PARAMETERS ----------------------- ###
### ---------------------------------------------------------------- ###

path2CLEO = sys.argv[1]  # absolute or relative path of CLEO directory
path2build = sys.argv[2] # absolute or relative path to build directory
isfigures = [True, True] # booleans for [making+showing, saving] figures

### --- essential paths and filenames --- ###
# where to find constants and config files
constsfile = path2CLEO+"/libs/claras_SDconstants.hpp"
configfile = path2CLEO+"/src/config/config.txt"

# where to write input files (note these should match configfile)
gridfile = path2build+"/share/dimlessGBxboundaries.dat"
thermofile =  path2build+"/share/dimlessthermodynamics.dat"
initSDsfile = path2build+"/share/dimlessSDsinit.dat"

# directory to save figures of initial conditions in
savefigpath = path2build+"/bin/"


### --- settings for gridbox boundaries --- ###
zgrid = [0, 1500, 50]     # evenly spaced zhalf coords [zmin, zmax, zdelta] [m]
xgrid = [0, 1500, 50]     # evenly spaced xhalf coords [m]
ygrid = np.array([0, 50])  # array of yhalf coords [m]


### --- settings for 2D Thermodyanmics --- ###
PRESS0 = 101500   # [Pa]
THETA = 289       # [K]
qvap = 0.0075     # [Kg/Kg]
qcond = 0.0       # [Kg/Kg]
WMAX = 0.6        # [m/s]
Zlength = 1500    # [m]
Xlength = 1500    # [m]
VVEL = None       # [m/s]


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
gen = gthermo.ConstHydrostaticAdiabat(configfile, constsfile, PRESS0, 
                                      THETA, qvap, qcond, WMAX, 
                                      Zlength, Xlength, VVEL)
cthermo.write_thermodynamics_binary(thermofile, gen, configfile,
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
  if isfigures[1]:
    Path(savefigpath).mkdir(exist_ok=True) 
  rgrid.plot_gridboxboundaries(constsfile, gridfile,
                               savefigpath, isfigures[1])
  rthermo.plot_thermodynamics(constsfile, configfile, gridfile,
                              thermofile, savefigpath, isfigures[1])
  rsupers.plot_initdistribs(configfile, constsfile, initSDsfile,
                            gridfile, savefigpath, isfigures[1]) 
### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###