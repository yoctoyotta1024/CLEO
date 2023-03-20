import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

from pySD.initsuperdropsbinary_src import *

### path and filenames
abspath = "/Users/yoctoyotta1024/Documents/b1_springsummer2023/CLEO/"
constsfile = abspath+"libs/claras_SDconstants.hpp"
configfile = abspath+"src/config/config.txt"

initSDsfilepath = abspath+"build/share/"
initSDsfile = initSDsfilepath+"dimlessSDsinit.dat"
binpath = abspath+"build/bin/"

### booleans for [making+saving, showing] figures
isfigures = [True, True]

### ------------ Choice of Superdroplet Radii Generator ------------ ###

monor                = 1e-6                        # all SDs have this same radius [m]
radiigen  = initattributes.MonoAttrsGen(monor)     # all SDs have the same dryradius [m]

rspan                = [1e-8, 9.1e-5]                # max and min range of radii to sample [m]
randomr              = False                        # sample radii range randomly or not
radiigen = initattributes.SampleDryradiiGen(rspan, randomr) # radii are sampled from rspan [m]


### ------ Choice of Droplet Radius Probability Distribution ------- ###

dirac0               = 1e-6                        # radius in sample closest to this value is dirac delta peak
numconc              = 1e9                         # total no. conc of real droplets [m^-3]
radiiprobdist = radiiprobdistribs.DiracDelta(dirac0)

# geomeans           = [0.075e-6]                  # lnnormal modes' geometric mean droplet radius [m] 
# geosigs            = [1.5]                       # lnnormal modes' geometric standard deviation
# scalefacs          = [1e9]                       # relative heights of modes         
geomeans             = [0.02e-6, 0.2e-6, 3.5e-6]               
geosigs              = [1.55, 2.3, 2]                    
scalefacs            = [1e6, 0.3e6, 0.025e6]   
numconc = np.sum(scalefacs) 
radiiprobdist = radiiprobdistribs.LnNormal(geomeans, geosigs, scalefacs)
 
volexpr0             = 30.531e-6                   # peak of volume exponential distribution [m]
numconc              = 2**(23)                     # total no. conc of real droplets [m^-3]
radiiprobdist = radiiprobdistribs.VolExponential(volexpr0, rspan)

### ---------- Choice of Superdroplet Coord3 Generator ------------- ###

#coord3gen            = None                        # do not generate superdroplet coord3s

# monocoord3           = 1000                        # all SDs have this same coord3 [m] 
# coord3gen = initattributes.MonoAttrsGen(monocoord3)

coord3span           = [0, 5000]                # max and min range of coord3 to sample [m]                 
randomcoord3         = True                        # sample coord3 range randomly or not
coord3gen = initattributes.SampleCoord3Gen(coord3span, randomcoord3)

### ---------------------------------------------------------------- ###


try:
  samplevol          = float(sys.argv[1])          # volume of droplet sample region [m^3]
except:
  errmsg = "please run script with volume of droplet sample region in m^3"
  raise ValueError(errmsg)

Path(binpath).mkdir(exist_ok=True) 
Path(initSDsfilepath).mkdir(exist_ok=True) 
initattrs = initattributes.InitAttributes(radiigen, radiiprobdist, coord3gen, numconc, samplevol)
create_initsuperdrops.write_initsuperdrops_binary(initSDsfile, initattrs, 
                                                  configfile, constsfile)

if isfigures[0]:
    read_initsuperdrops.plot_initdistribs(configfile, constsfile, initSDsfile,
                                          samplevol, binpath, isfigures[1])