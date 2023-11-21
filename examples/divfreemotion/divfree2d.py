'''
----- CLEO -----
File: divfree2d.py
Project: divfreemotion
Created Date: Friday 17th November 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
Last Modified: Tuesday 21st November 2023
Modified By: CB
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
Copyright (c) 2023 MPI-M, Clara Bayley
-----
File Description:
Script compiles and runs CLEO adia0D to create the
data and plots as in Shima et al. 2009 to show
comparision of SDM 0D box model of collision-
coalescence with Golovin's analytical solution
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

from pySD.sdmout_src import *
from pySD.initsuperdropsbinary_src import *
from pySD.gbxboundariesbinary_src import read_gbxboundaries as rgrid
from pySD.gbxboundariesbinary_src import create_gbxboundaries as cgrid

### ---------------------------------------------------------------- ###
### ----------------------- INPUT PARAMETERS ----------------------- ###
### ---------------------------------------------------------------- ###

### --- essential paths and filenames --- ###
# path and filenames for creating initial SD conditions
constsfile = path2CLEO+"/libs/cleoconstants.hpp"
binpath = path2build+"/bin/"
sharepath = path2build+"/share/"
gridfile = sharepath+"divfree2d_dimlessGBxboundaries.dat"
initSDsfile = sharepath+"divfree2d_dimlessSDsinit.dat"
thermofile =  sharepath+"/divfree2d_dimlessthermo.dat"

# path and file names for plotting results
setupfile = binpath+"divfree2d_setup.txt"
dataset = binpath+"divfree2d_sol.zarr"

# booleans for [making, showing] initialisation figures
isfigures = [True, False]

### --- settings for 2-D gridbox boundaries --- ###
zgrid = [0, 1500, 100]     # evenly spaced zhalf coords [zmin, zmax, zdelta] [m]
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
  Path(sharepath).mkdir(exist_ok=True) 
  Path(binpath).mkdir(exist_ok=True) 

### 1. create files for gridbox boundaries and initial SD conditions
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

# 2. compile and the run model
os.chdir(path2build)
os.system('pwd')
os.system('rm -rf '+dataset)
os.system('make -j 64 divfree2D')
# os.system('make clean && make -j 64 gol0D')
executable = path2build+'/examples/divfreemotion/src/divfree2D'
os.system(executable + ' ' + configfile)

# 3. load and plot results
# read in constants and intial setup from setup .txt file
config = pysetuptxt.get_config(setupfile, nattrs=3, isprint=True)
consts = pysetuptxt.get_consts(setupfile, isprint=True)
gbxs = pygbxsdat.get_gridboxes(gridfile, consts["COORD0"], isprint=True)

time = pyzarr.get_time(dataset).secs
sddata = pyzarr.get_supers(dataset, consts)

# 4. plot results
savename = binpath + "divfree2d_validation.png"