"""
----- CLEO -----
File: cuspbifurc.py
Project: adiabaticparcel
Created Date: Friday 17th November 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
Last Modified: Friday 3rd May 2024
Modified By: CB
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
Copyright (c) 2023 MPI-M, Clara Bayley
-----
File Description:
Script generate input files, runs CLEO adia0d executable to
create data and then creates plots for adiabatic parcel example
similar to Figure 5 of "On the CCN (de)activation nonlinearities"
S. Arabas and S. Shima 2017 to show example of cusp birfucation for
0D adiabatic parcel expansion and contraction.
Note: SD(M) = superdroplet (model)
"""

import os
import shutil
import subprocess
import sys
import numpy as np
from pathlib import Path

path2CLEO = Path(sys.argv[1])
path2build = Path(sys.argv[2])
config_filename = Path(sys.argv[3])

sys.path.append(str(path2CLEO))  # imports from pySD
sys.path.append(
    str(path2CLEO / "examples" / "exampleplotting")
)  # imports from example plots package


from plotssrc import pltsds, as2017fig
from pySD import geninitconds
from pySD.sdmout_src import pyzarr, pysetuptxt, pygbxsdat, sdtracing
from pySD.initsuperdropsbinary_src import rgens, dryrgens, probdists, attrsgen


############### INPUTS ##################
# path and filenames for creating SD initial conditions and for running model
constants_filename = path2CLEO / "libs" / "cleoconstants.hpp"
binpath = path2build / "bin"
sharepath = path2build / "share"
initsupers_filename = sharepath / "cuspbifurc_dimlessSDsinit.dat"
grid_filename = sharepath / "cuspbifurc_dimlessGBxboundaries.dat"

# path and file names for plotting results
setupfile = binpath / "cuspbifurc_setup.txt"
dataset = binpath / "cuspbifurc_sol.zarr"

# booleans for [making, saving] initialisation figures
isfigures = [True, True]
savefigpath = binpath

# settings for 0D Model (number of SD and grid coordinates)
nsupers = {0: 1}
coord_params = ["false"]
zgrid = np.asarray([0, 100])
xgrid = np.asarray([0, 100])
ygrid = np.asarray([0, 100])

# settings for monodisperse droplet radii
# numconc = total no. concentration of droplets [m^-3]
numconc = 0.5e9
# monor = dry radius of all droplets [m]
monor = 0.025e-6

# monodisperse droplet radii probability distribution
radiigen = rgens.MonoAttrGen(monor)
dryradiigen = dryrgens.ScaledRadiiGen(1.0)
xiprobdist = probdists.DiracDelta(monor)

# do not generate SD coords
coord3gen = None
coord1gen = None
coord2gen = None


def displacement(time, w_avg, thalf):
    """displacement z given velocity, w, is sinusoidal
    profile: w = w_avg * pi/2 * np.sin(np.pi * t/thalf)
    where wmax = pi/2*w_avg and tauhalf = thalf/pi."""

    zmax = w_avg / 2 * thalf
    z = zmax * (1 - np.cos(np.pi * time / thalf))
    return z


############### RUN EXAMPLE ##################
### ensure build, share and bin directories exist
if path2CLEO == path2build:
    raise ValueError("build directory cannot be CLEO")
else:
    path2build.mkdir(exist_ok=True)
    sharepath.mkdir(exist_ok=True)
    binpath.mkdir(exist_ok=True)
    if isfigures[1]:
        savefigpath.mkdir(exist_ok=True)

###  delete any exisitng initial conditions
shutil.rmtree(grid_filename, ignore_errors=True)
shutil.rmtree(initsupers_filename, ignore_errors=True)

### create files (and plots) for gridbox boundaries and initial SD conditions
geninitconds.generate_gridbox_boundaries(
    grid_filename,
    zgrid,
    xgrid,
    ygrid,
    constants_filename,
    isprintinfo=True,
    isfigures=isfigures,
    savefigpath=savefigpath,
)

initattrsgen = attrsgen.AttrsGenerator(
    radiigen, dryradiigen, xiprobdist, coord3gen, coord1gen, coord2gen
)
geninitconds.generate_initial_superdroplet_conditions(
    initattrsgen,
    initsupers_filename,
    config_filename,
    constants_filename,
    grid_filename,
    nsupers,
    numconc,
    isprintinfo=True,
    isfigures=isfigures,
    savefigpath=savefigpath,
    gbxs2plt="all",
)

### run model
os.chdir(path2build)
subprocess.run(["pwd"])
shutil.rmtree(dataset, ignore_errors=True)  # delete any existing dataset
executable = path2build / "examples" / "adiabaticparcel" / "src" / "adia0d"
print("Executable: " + str(executable))
print("Config file: " + str(config_filename))
subprocess.run([executable, config_filename])

### load results
# read in constants and intial setup from setup .txt file
config = pysetuptxt.get_config(setupfile, nattrs=3, isprint=True)
consts = pysetuptxt.get_consts(setupfile, isprint=True)
gbxs = pygbxsdat.get_gridboxes(grid_filename, consts["COORD0"], isprint=True)

# read in output Xarray data
thermo = pyzarr.get_thermodata(dataset, config["ntime"], gbxs["ndims"], consts)
supersat = thermo.supersaturation()
time = pyzarr.get_time(dataset).secs
sddata = pyzarr.get_supers(dataset, consts)
zprof = displacement(time, config["W_avg"], config["TAU_half"])

### plot results
# sample drops to plot from whole range of SD ids
sample = [0, int(config["maxnsupers"])]
radii = sdtracing.attribute_for_superdroplets_sample(
    sddata, "radius", minid=sample[0], maxid=sample[1]
)
savename = savefigpath / "cuspbifurc_SDgrowth.png"
pltsds.individ_radiusgrowths_figure(time, radii, savename=savename)

attrs = ["radius", "xi", "msol"]
sd0 = sdtracing.attributes_for1superdroplet(sddata, 0, attrs)
numconc = np.sum(sddata["xi"][0]) / gbxs["domainvol"] / 1e6  # [/cm^3]

savename2 = savefigpath / "cuspbifurc_validation.png"
as2017fig.arabas_shima_2017_fig(
    time,
    zprof,
    sd0["radius"],
    sd0["msol"],
    thermo.temp[:, 0, 0, 0],
    supersat[:, 0, 0, 0],
    sddata.IONIC,
    sddata.MR_SOL,
    config["W_avg"],
    numconc,
    savename2,
)
