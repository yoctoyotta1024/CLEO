import sys
import numpy as np
from pathlib import Path

from pySD.initsuperdropsbinary_src import initattributes as iSDs
from pySD.initsuperdropsbinary_src import radiiprobdistribs as rprobs
from pySD.initsuperdropsbinary_src import create_initsuperdrops as csupers 
from pySD.initsuperdropsbinary_src import read_initsuperdrops as rsupers 

### ----------------------- INPUT PARAMETERS ----------------------- ###
### --- absolute or relative paths for --- ###
### ---   build and CLEO directories --- ###
path2CLEO = sys.argv[1]
path2build = sys.argv[2]
configfile = sys.argv[3]

### booleans for [making+showing, saving] figures
isfigures = [True, True]
gbxs2plt = 0 # indexes of GBx index of SDs to plot (nb. "all" can be very slow)

### essential paths and filenames
constsfile = path2CLEO+"libs/claras_SDconstants.hpp"
binariespath = path2build+"/share/"
savefigpath = path2build+"/bin/"

gridfile =  binariespath+"/dimlessGBxboundaries.dat" # note this should match config.txt
initSDsfile = binariespath+"/dimlessSDsinit.dat" # note this should match config.txt

### --- Number of Superdroplets per Gridbox --- ###
### ---        (an int or dict of ints)     --- ###
# zlim = 1000
# npergbx = 8
# nsupers = iSDs.nsupers_at_domain_base(gridfile, constsfile, npergbx, zlim)
nsupers = 1024
### ------------------------------------------- ###

### --- Choice of Superdroplet Radii Generator --- ###
# monor                = 0.04910258806                  # all SDs have this same radius [m]
# radiigen  = iSDs.MonoAttrsGen(monor)                  # all SDs have the same dryradius [m]

rspan                = [1e-8, 9e-5]                 # min and max range of radii to sample [m]
randomr              = True                         # sample radii range randomly or not
radiigen = iSDs.SampleDryradiiGen(rspan, randomr)   # radii are sampled from rspan [m]
### ---------------------------------------------- ###

### --- Choice of Droplet Radius Probability Distribution --- ###
# dirac0               = monor               # radius in sample closest to this value is dirac delta peak
# numconc              = 1e6                         # total no. conc of real droplets [m^-3]
# radiiprobdist = rprobs.DiracDelta(dirac0)

# geomeans           = [0.075e-6]                  # lnnormal modes' geometric mean droplet radius [m] 
# geosigs            = [1.5]                       # lnnormal modes' geometric standard deviation
# scalefacs          = [1e9]                       # relative heights of modes         
# geomeans             = [0.02e-6, 0.2e-6, 3.5e-6]               
# geosigs              = [1.55, 2.3, 2]                    
# scalefacs            = [1e6, 0.3e6, 0.025e6]   
# geomeans             = [0.02e-6, 0.15e-6]               
# geosigs              = [1.4, 1.6]                    
# scalefacs            = [6e6, 4e6]   
# numconc = np.sum(scalefacs)
# radiiprobdist = rprobs.LnNormal(geomeans, geosigs, scalefacs)
 
volexpr0             = 30.531e-6                   # peak of volume exponential distribution [m]
numconc              = 2**(23)                     # total no. conc of real droplets [m^-3]
radiiprobdist = rprobs.VolExponential(volexpr0, rspan)

# reff                 = 7e-6                     # effective radius [m]
# nueff                = 0.08                     # effective variance 
# # radiiprobdist = rprobs.ClouddropsHansenGamma(reff, nueff)
# rdist1 = rprobs.ClouddropsHansenGamma(reff, nueff)
# nrain                = 3000                         # raindrop concentration [m^-3]
# qrain                = 0.9                          # rainwater content [g/m^3]
# dvol                 = 8e-4                         # mean volume diameter [m]
# # radiiprobdist = rprobs.RaindropsGeoffroyGamma(nrain, qrain, dvol)
# rdist2 = rprobs.RaindropsGeoffroyGamma(nrain, qrain, dvol)
# numconc = 1e9 # [m^3]
# distribs = [rdist1, rdist2]
# scalefacs = [1000, 1]
# radiiprobdist = rprobs.CombinedRadiiProbDistribs(distribs, scalefacs)

### --------------------------------------------------------- ###

### --- Choice of Superdroplet Coord3 Generator --- ###
# monocoord3           = 1000                        # all SDs have this same coord3 [m] 
# coord3gen            = iSDs.MonoCoordGen(monocoord3)
# coord3gen            = iSDs.SampleCoordGen(True) # sample coord3 range randomly or not
coord3gen            = None                        # do not generate superdroplet coord3s
### ----------------------------------------------- ###

### --- Choice of Superdroplet Coord1 Generator --- ###
# monocoord1           = 200                        # all SDs have this same coord1 [m] 
# coord1gen            = iSDs.MonoCoordGen(monocoord1)
# coord1gen            = iSDs.SampleCoordGen(True) # sample coord1 range randomly or not
coord1gen            = None                        # do not generate superdroplet coord1s
### ----------------------------------------------- ###

### --- Choice of Superdroplet Coord2 Generator --- ###
# monocoord2           = 1000                        # all SDs have this same coord2 [m] 
# coord2gen            = iSDs.MonoCoordGen(monocoord2)
# coord2gen            = iSDs.SampleCoordGen(True) # sample coord1 range randomly or not
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
initattrsgen = iSDs.InitManyAttrsGen(radiigen, radiiprobdist,
                                      coord3gen, coord1gen, coord2gen)
csupers.write_initsuperdrops_binary(initSDsfile, initattrsgen, 
                                    configfile, constsfile,
                                    gridfile, nsupers, numconc)

### plot initial superdrops binary
if isfigures[0]:
    if isfigures[1]:
        Path(savefigpath).mkdir(exist_ok=True) 
    rsupers.plot_initGBxsdistribs(configfile, constsfile, initSDsfile,
                                   gridfile, savefigpath, isfigures[1],
                                   gbxs2plt)
### ---------------------------------------------------------------- ###