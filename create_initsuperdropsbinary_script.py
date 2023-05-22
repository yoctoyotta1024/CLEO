import sys
import numpy as np
from pathlib import Path

import pySD.initsuperdropsbinary_src as iSDs
from pySD.gbxboundariesbinary_src.read_gbxboundaries import get_domainvol_from_gridfile

### ----------------------- INPUT PARAMETERS ----------------------- ###
### --- absolute or relative paths for --- ###
### ---   build and CLEO directories --- ###
path2CLEO = sys.argv[1]
path2build = sys.argv[2]

### booleans for [making+showing, saving] figures
# isfigures = [True, True]
isfigures = [False, False]

### essential paths and filenames
constsfile = path2CLEO+"libs/claras_SDconstants.hpp"
configfile = path2CLEO+"src/config/config.txt"
binariespath = path2build+"/share/"
savefigpath = path2build+"/bin/"

gridfile =  binariespath+"/dimlessGBxboundaries.dat" # note this should match config.txt
initSDsfile = binariespath+"/dimlessSDsinit.dat" # note this should match config.txt

### --- Number of Superdroplets per Gridbox --- ###
### ---        (an int or dict of ints)     --- ###
zlim = 750
npergbx = 1024
nsupers = iSDs.initattributes.nsupers_at_domain_base(gridfile, constsfile, npergbx, zlim)
# nsupers = 1024
### ------------------------------------------- ###

### --- Choice of Superdroplet Radii Generator --- ###
# monor                = 1e-6                        # all SDs have this same radius [m]
# radiigen  = iSDs.initattributes.MonoAttrsGen(monor)     # all SDs have the same dryradius [m]

# rspan                = [1e-8, 9e-5]                # min and max range of radii to sample [m]
rspan                = [3e-9, 3e-6]                # min and max range of radii to sample [m]
randomr              = True                        # sample radii range randomly or not
radiigen = iSDs.initattributes.SampleDryradiiGen(rspan, randomr) # radii are sampled from rspan [m]
### ---------------------------------------------- ###

### --- Choice of Droplet Radius Probability Distribution --- ###
# dirac0               = 1e-6                        # radius in sample closest to this value is dirac delta peak
# numconc              = 1e9                         # total no. conc of real droplets [m^-3]
# radiiprobdist = iSDs.radiiprobdistribs.DiracDelta(dirac0)

# geomeans           = [0.075e-6]                  # lnnormal modes' geometric mean droplet radius [m] 
# geosigs            = [1.5]                       # lnnormal modes' geometric standard deviation
# scalefacs          = [1e9]                       # relative heights of modes         
# geomeans             = [0.02e-6, 0.2e-6, 3.5e-6]               
# geosigs              = [1.55, 2.3, 2]                    
# scalefacs            = [1e6, 0.3e6, 0.025e6]   
geomeans             = [0.02e-6, 0.15e-6]               
geosigs              = [1.4, 1.6]                    
scalefacs            = [6e6, 4e6]   
numconc = np.sum(scalefacs)
radiiprobdist = iSDs.radiiprobdistribs.LnNormal(geomeans, geosigs, scalefacs)
 
# volexpr0             = 30.531e-6                   # peak of volume exponential distribution [m]
# numconc              = 2**(23)                     # total no. conc of real droplets [m^-3]
# radiiprobdist = iSDs.radiiprobdistribs.VolExponential(volexpr0, rspan)
### --------------------------------------------------------- ###

### --- Choice of Superdroplet Coord3 Generator --- ###
# monocoord3           = 1000                        # all SDs have this same coord3 [m] 
# coord3gen            = iSDs.initattributes.MonoCoordGen(monocoord3)
coord3gen            = iSDs.initattributes.SampleCoordGen(True) # sample coord3 range randomly or not
# coord3gen            = None                        # do not generate superdroplet coord3s
### ----------------------------------------------- ###

### --- Choice of Superdroplet Coord1 Generator --- ###
# monocoord1           = 200                        # all SDs have this same coord1 [m] 
# coord1gen            = iSDs.initattributes.MonoCoordGen(monocoord1)
coord1gen            = iSDs.initattributes.SampleCoordGen(True) # sample coord1 range randomly or not
# coord1gen            = None                        # do not generate superdroplet coord1s
### ----------------------------------------------- ###

### --- Choice of Superdroplet Coord2 Generator --- ###
# monocoord2           = 1000                        # all SDs have this same coord2 [m] 
# coord2gen            = iSDs.initattributes.MonoCoordGen(monocoord2)
# coord2gen            = iSDs.initattributes.SampleCoordGen(True) # sample coord1 range randomly or not
coord2gen            = None                        # do not generate superdroplet coord2s
### ----------------------------------------------- ###

### ---------------------------------------------------------------- ###


### -------------------- BINARY FILE GENERATION--------------------- ###
### ensure build, share and bin directories exist
if path2CLEO == path2build:
  raise ValueError("build directory cannot be CLEO")
else:
  Path(path2build).mkdir(exist_ok=True) 
  Path(binariespath).mkdir(exist_ok=True) 

### write initial superdrops binary
initattrsgen = iSDs.initattributes.InitManyAttrsGen(radiigen, radiiprobdist,
                                               coord3gen, coord1gen, coord2gen)
iSDs.create_initsuperdrops.write_initsuperdrops_binary(initSDsfile, initattrsgen, 
                                                  configfile, constsfile,
                                                  gridfile, nsupers, numconc)

### plot initial superdrops binary
if isfigures[0]:
    if isfigures[1]:
        Path(savefigpath).mkdir(exist_ok=True) 
    iSDs.read_initsuperdrops.plot_initdistribs(configfile, constsfile, initSDsfile,
                                          gridfile, savefigpath, isfigures[1])
### ---------------------------------------------------------------- ###