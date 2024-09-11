"""
----- CLEO -----
File: rainshaft1d.py
Project: rainshaft1d
Created Date: Friday 17th November 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
Last Modified: Friday 12th April 2024
Modified By: CB
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
Copyright (c) 2023 MPI-M, Clara Bayley
-----
File Description:
Script generates input files, then runs CLEO executable "rshaft1d" to create the
data which is then plotted to demonstrate precipitation example in 1-D rainshaft
with constant thermodynamics read from a file.
"""

import os
import sys
import numpy as np
import random
from pathlib import Path

path2CLEO = sys.argv[1]
path2build = sys.argv[2]
configfile = sys.argv[3]

sys.path.append(path2CLEO)  # for imports from pySD package
sys.path.append(
    path2CLEO + "/examples/exampleplotting/"
)  # for imports from example plotting package

from plotssrc import pltsds, pltmoms, animations
from pySD.sdmout_src import pyzarr, pysetuptxt, pygbxsdat
from pySD.gbxboundariesbinary_src import read_gbxboundaries as rgrid
from pySD.gbxboundariesbinary_src import create_gbxboundaries as cgrid
from pySD.initsuperdropsbinary_src import crdgens, rgens, dryrgens, probdists, attrsgen
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
constsfile = path2CLEO + "/libs/cleoconstants.hpp"
binpath = path2build + "/bin/"
sharepath = path2build + "/share/"
gridfile = sharepath + "rain1d_dimlessGBxboundaries.dat"
initSDsfile = sharepath + "rain1d_dimlessSDsinit.dat"
thermofile = sharepath + "rain1d_dimlessthermo.dat"

# path and file names for plotting results
setupfile = binpath + "rain1d_setup.txt"
dataset = binpath + "rain1d_sol.zarr"

### --- plotting initialisation figures --- ###
isfigures = [True, True]  # booleans for [making, saving] initialisation figures
savefigpath = path2build + "/bin/"  # directory for saving figures
SDgbxs2plt = list(range(39, 124))
SDgbxs2plt = [random.choice(SDgbxs2plt)]  # choose random gbx from list to plot

### --- settings for 1-D gridbox boundaries --- ###
zgrid = [0, 2500, 20]  # evenly spaced zhalf coords [zmin, zmax, zdelta] [m]
xgrid = np.array([0, 20])  # array of xhalf coords [m]
ygrid = np.array([0, 20])  # array of yhalf coords [m]

### --- settings for 1-D Thermodynamics --- ###
PRESS0 = 101315  # [Pa]
TEMP0 = 297.9  # [K]
qvap0 = 0.016  # [Kg/Kg]
Zbase = 800  # [m]
TEMPlapses = [9.8, 6.5]  # -dT/dz [K/km]
qvaplapses = [2.97, "saturated"]  # -dvap/dz [g/Kg km^-1]
qcond = 0.0  # [Kg/Kg]
WVEL = 4.0  # [m/s]
Wlength = (
    1000  # [m] use constant W (Wlength=0.0), or sinusoidal 1-D profile below cloud base
)

### --- settings for initial superdroplets --- ###
# initial superdroplet coordinates
zlim = 800  # min z coord of superdroplets [m]
npergbx = 256  # number of superdroplets per gridbox

# initial superdroplet radii (and implicitly solute masses)
rspan = [3e-9, 5e-5]  # min and max range of radii to sample [m]
dryr_sf = 1.0  # dryradii are 1/sf of radii [m]

# settings for initial superdroplet multiplicies
geomeans = [0.02e-6, 0.2e-6, 3.5e-6]
geosigs = [1.55, 2.3, 2]
scalefacs = [1e6, 0.3e6, 0.025e6]
numconc = np.sum(scalefacs) * 1000
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

### --- delete any existing initial conditions --- ###
os.system("rm " + gridfile)
os.system("rm " + initSDsfile)
os.system("rm " + thermofile[:-4] + "*")

### ----- write gridbox boundaries binary ----- ###
cgrid.write_gridboxboundaries_binary(gridfile, zgrid, xgrid, ygrid, constsfile)
rgrid.print_domain_info(constsfile, gridfile)

### ----- write thermodynamics binaries ----- ###
thermodyngen = thermogen.ConstHydrostaticLapseRates(
    configfile,
    constsfile,
    PRESS0,
    TEMP0,
    qvap0,
    Zbase,
    TEMPlapses,
    qvaplapses,
    qcond,
    WVEL,
    None,
    None,
    Wlength,
)
cthermo.write_thermodynamics_binary(
    thermofile, thermodyngen, configfile, constsfile, gridfile
)

### ----- write initial superdroplets binary ----- ###
nsupers = crdgens.nsupers_at_domain_top(gridfile, constsfile, npergbx, zlim)
coord3gen = crdgens.SampleCoordGen(True)  # sample coord3 randomly
coord1gen = None  # do not generate superdroplet coord2s
coord2gen = None  # do not generate superdroplet coord2s

xiprobdist = probdists.LnNormal(geomeans, geosigs, scalefacs)
radiigen = rgens.SampleLog10RadiiGen(rspan)
dryradiigen = dryrgens.ScaledRadiiGen(dryr_sf)

initattrsgen = attrsgen.AttrsGenerator(
    radiigen, dryradiigen, xiprobdist, coord3gen, coord1gen, coord2gen
)
csupers.write_initsuperdrops_binary(
    initSDsfile, initattrsgen, configfile, constsfile, gridfile, nsupers, numconc
)

### ----- show (and save) plots of binary file data ----- ###
if isfigures[0]:
    if isfigures[1]:
        Path(savefigpath).mkdir(exist_ok=True)
    rgrid.plot_gridboxboundaries(constsfile, gridfile, savefigpath, isfigures[1])
    rthermo.plot_thermodynamics(
        constsfile, configfile, gridfile, thermofile, savefigpath, isfigures[1]
    )
    rsupers.plot_initGBxs_distribs(
        configfile,
        constsfile,
        initSDsfile,
        gridfile,
        savefigpath,
        isfigures[1],
        SDgbxs2plt,
    )
### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###

### ---------------------------------------------------------------- ###
### ---------------------- RUN CLEO EXECUTABLE --------------------- ###
### ---------------------------------------------------------------- ###
os.chdir(path2build)
os.system("pwd")
os.system("rm -rf " + dataset)  # delete any existing dataset
executable = path2build + "/examples/rainshaft1d/src/rshaft1d"
print("Executable: " + executable)
print("Config file: " + configfile)
os.system(executable + " " + configfile)
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
pltsds.plot_randomsample_superdrops(
    time, sddata, config["maxnsupers"], nsample, savename=savename
)

### ----- plot 1-D .gif animations ----- ###
nframes = len(time.mins)
mom2ani = np.sum(massmoms.nsupers, axis=(1, 2))
xlims = [0, np.amax(mom2ani)]
xlabel = "number of super-droplets"
savename = savefigpath + "rain1d_nsupers1d"
animations.animate1dprofile(
    gbxs,
    mom2ani,
    time.mins,
    nframes,
    xlabel=xlabel,
    xlims=xlims,
    color="green",
    saveani=True,
    savename=savename,
    fps=5,
)

nframes = len(time.mins)
norm = gbxs["gbxvols"] * 1e6  # volume [cm^3]
mom2ani = np.sum(massmoms.mom0 / norm[None, :], axis=(1, 2))
xlims = [0, np.amax(mom2ani)]
xlabel = "number concentration /cm$^{-3}$"
savename = savefigpath + "rain1d_numconc1d"
animations.animate1dprofile(
    gbxs,
    mom2ani,
    time.mins,
    nframes,
    xlabel=xlabel,
    xlims=xlims,
    color="green",
    saveani=True,
    savename=savename,
    fps=5,
)

nframes = len(time.mins)
norm = gbxs["gbxvols"]  # volume [m^3]
mom2ani = np.sum(massmoms.mom1 / norm[None, :], axis=(1, 2))
xlims = [0, np.amax(mom2ani)]
xlabel = "mass concentration /g m$^{-3}$"
savename = savefigpath + "rain1d_massconc1d"
animations.animate1dprofile(
    gbxs,
    mom2ani,
    time.mins,
    nframes,
    xlabel=xlabel,
    xlims=xlims,
    color="green",
    saveani=True,
    savename=savename,
    fps=5,
)
### ------------------------------------------------------------ ###
### ------------------------------------------------------------ ###
