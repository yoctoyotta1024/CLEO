'''
----- CLEO -----
File: yac1_fromfile.py
Project: yac1_fromfile
Created Date: Friday 17th November 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
Last Modified: Friday 22nd March 2024
Modified By: CB
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
Copyright (c) 2023 MPI-M, Clara Bayley
-----
File Description:
Script compiles and runs CLEO for 3D example with time varying thermodynamics
read from binary files to test that YAC can send the data to CLEO correctly.
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
# for imports from example plotting package
sys.path.append(path2CLEO+"/examples/exampleplotting/")

from plotssrc import pltsds, pltmoms
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
gridfile = sharepath+"yac1_dimlessGBxboundaries.dat"
initSDsfile = sharepath+"yac1_dimlessSDsinit.dat"
thermofile = sharepath+"/yac1_dimlessthermo.dat"

# path and file names for plotting results
setupfile = binpath+"yac1_setup.txt"
dataset = binpath+"yac1_sol.zarr"

### --- plotting initialisation figures --- ###
# booleans for [making, saving] initialisation figures
isfigures = [True, True]
savefigpath = path2build+"/bin/"  # directory for saving figures
SDgbxs2plt = [0]  # gbxindex of SDs to plot (nb. "all" can be very slow)

### --- settings for 2-D gridbox boundaries --- ###
zgrid = [0, 1500, 60]           # evenly spaced zhalf coords [zmin, zmax, zdelta] [m]
xgrid = [0, 1500, 60]           # evenly spaced xhalf coords [m]
ygrid = np.array([0, 10, 20])   # array of yhalf coords [m]

### --- settings for initial superdroplets --- ###
# settings for initial superdroplet coordinates
zlim = 1000             # max z coord of superdroplets
npergbx = 2             # number of superdroplets per gridbox

monor = 1e-6            # all SDs have this same radius [m]
dryr_sf = 1.0           # scale factor for dry radii [m]
numconc = 5e8           # total no. conc of real droplets [m^-3]
randcoord = False       # sample SD spatial coordinates randomly or not

### --- settings for 2D Thermodynamics --- ###
PRESS0 = 100000 # [Pa]
THETA = 298.15  # [K]
qcond = 0.0     # [Kg/Kg]
WMAX = 0.6      # [m/s]
VVEL = 2.0      # [m/s]
Zlength = 1500  # [m]
Xlength = 1500  # [m]
qvapmethod = "sratio"
Zbase = 750     # [m]
moistlayer = False
sratios = [1.0, 1.0]  # s_ratio [below, above] Zbase
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

### ----- write thermodynamics binaries ----- ###
thermodyngen = thermogen.ConstDryHydrostaticAdiabat(configfile, constsfile, PRESS0,
                                        THETA, qvapmethod, sratios, Zbase,
                                        qcond, WMAX, Zlength, Xlength,
                                        VVEL, moistlayer)
cthermo.write_thermodynamics_binary(thermofile, thermodyngen, configfile,
                                    constsfile, gridfile)

### ----- write initial superdroplets binary ----- ###
nsupers = crdgens.nsupers_at_domain_base(gridfile, constsfile, npergbx, zlim)
radiigen  =  rgens.MonoAttrGen(monor)                 # all SDs have the same radius [m]
dryradiigen =  dryrgens.ScaledRadiiGen(dryr_sf)       # dryradii are 1/sf of radii [m]
coord3gen =  crdgens.SampleCoordGen(randcoord)        # (not) random coord3 of SDs
coord1gen =  crdgens.SampleCoordGen(randcoord)        # (not) random coord1 of SDs
coord2gen =  crdgens.SampleCoordGen(randcoord)        # (not) random coord2 of SDs
xiprobdist = probdists.DiracDelta(monor)              # monodisperse droplet probability distrib

initattrsgen = attrsgen.AttrsGenerator(radiigen, dryradiigen, xiprobdist,
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
    rsupers.plot_initGBxs_distribs(configfile, constsfile, initSDsfile,
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
os.system('make clean && make -j 64 yac1')
executable = path2build+'/examples/yac_examples/yac1_fromfile/src/yac1'
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
totnsupers = pyzarr.get_totnsupers(dataset)

# 4. plot results
savename = savefigpath + "yac1_totnsupers_validation.png"
pltmoms.plot_totnsupers(time, totnsupers, savename=savename)

nsample = 500
savename = savefigpath + "yac1_motion2d_validation.png"
pltsds.plot_randomsample_superdrops_2dmotion(sddata,
                                             config["totnsupers"],
                                             nsample,
                                             savename=savename,
                                             arrows=False)
### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###
