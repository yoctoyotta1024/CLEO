'''
----- CLEO -----
File: adiabaticparcel.py
Project: adiabaticparcel
Created Date: Friday 17th November 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
Last Modified: Friday 17th November 2023
Modified By: CB
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
Copyright (c) 2023 MPI-M, Clara Bayley
-----
File Description:
Script compiles and runs CLEO adia0D to
create data and plots similar to Figure 5 of
"On the CCN (de)activation nonlinearities"
S. Arabas and S. Shima 2017 to show
example of cusp birfucation for
0D adaibatic parcel exapansion and contraction.
Note: SD(M) = superdroplet (model)
'''

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

path2CLEO = sys.argv[1]
path2build = sys.argv[2]
configfile = sys.argv[3]

# for imports from pySD package
sys.path.append(path2CLEO)
from pySD.gbxboundariesbinary_src import create_gbxboundaries as cgrid
from pySD.gbxboundariesbinary_src import read_gbxboundaries as rgrid
from pySD.initsuperdropsbinary_src import *

############### INPUTS ##################
# path and filenames for creating SD initial conditions and for running model
constsfile = path2CLEO+"libs/cleoconstants.hpp"
binpath = path2build+"/bin/"
sharepath = path2build+"/share/"
initSDsfile = sharepath+"adai0D_dimlessSDsinit.dat"
gridfile = sharepath+"adia0D_dimlessGBxboundaries.dat"

# path and file names for plotting results
setupfile = binpath+"adia0Dsetup.txt"
dataset = binpath+"adia0Dsol.zarr"

# booleans for [making, showing] initialisation figures
isfigures = [True, True]

# settings for 0D Model (number of SD and grid coordinates)
nsupers = {0: 1}
coord_params = ["false"]
zgrid = np.asarray([0, 100])
xgrid = np.asarray([0, 100])
ygrid = np.asarray([0, 100])

# settings for monodisperse droplet radii
numconc = 0.5e9  # numconc = total no. concentration of droplets [m^-3]
monor = 0.025e-6  # monor = dry radius of all droplets [m]

# monodisperse droplet radii probability distribution
radiigen = initattributes.MonoAttrsGen(monor)
radiiprobdist = radiiprobdistribs.DiracDelta(monor)

# volume SD sample occupies (entire domain) [m^3]
samplevol = rgrid.calc_domainvol(zgrid, xgrid, ygrid)
coord3gen = None                        # do not generate SD coords
coord1gen = None
coord2gen = None

def displacement(time, w_avg, thalf):
  '''displacement z given velocity, w, is sinusoidal 
  profile: w = w_avg * pi/2 * np.sin(np.pi * t/thalf)
  where wmax = pi/2*w_avg and tauhalf = thalf/pi.'''

  zmax = w_avg / 2 * thalf
  z = zmax * (1 - np.cos(np.pi * time / thalf))
  return z

### 1. create files for gridbox boundaries and initial SD conditions
Path(binpath).mkdir(parents=True, exist_ok=True)
os.system("rm "+gridfile)
os.system("rm "+initSDsfile)
cgrid.write_gridboxboundaries_binary(gridfile, zgrid, xgrid,
                                     ygrid, constsfile)
rgrid.print_domain_info(constsfile, gridfile)

initattrsgen = initattributes.InitManyAttrsGen(radiigen, radiiprobdist,
                                               coord3gen, coord1gen, coord2gen)
create_initsuperdrops.write_initsuperdrops_binary(initSDsfile, initattrsgen,
                                                  configfile, constsfile,
                                                  gridfile, nsupers, numconc)
read_initsuperdrops.print_initSDs_infos(initSDsfile, configfile, constsfile, gridfile)

if isfigures[0]:
    rgrid.plot_gridboxboundaries(constsfile, gridfile,
                                 binpath, isfigures[1])
    read_initsuperdrops.plot_initGBxsdistribs(configfile, constsfile, initSDsfile,
                                              gridfile, binpath, isfigures[1], "all")
plt.close()