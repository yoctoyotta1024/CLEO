"""
Copyright (c) 2024 MPI-M, Clara Bayley


----- CLEO -----
File: shima2009.py
Project: boxmodelcollisions
Created Date: Friday 17th November 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
Last Modified: Sunday 16th June 2024
Modified By: CB
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
Script generates input files, runs CLEO 0-D box model executables for
collisions with selected collision kernels (e.g. Golovin's or Long's) to create data.
Then plots results comparable to Shima et al. 2009 Fig. 2
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
kernels = sys.argv[4:]

sys.path.append(str(path2CLEO))  # imports from pySD
sys.path.append(
    str(path2CLEO / "examples" / "exampleplotting")
)  # imports from example plots package


import attrgens_shima2009
from plotssrc import shima2009fig
from pySD import editconfigfile, geninitconds
from pySD.sdmout_src import pyzarr, pysetuptxt, pygbxsdat
from pySD.initsuperdropsbinary_src import rgens, attrsgen
from pySD.gbxboundariesbinary_src import read_gbxboundaries as rgrid

### ---------------------------------------------------------------- ###
### ----------------------- INPUT PARAMETERS ----------------------- ###
### ---------------------------------------------------------------- ###
### --- essential paths and filenames --- ###
# path and filenames for creating initial SD conditions
constants_filename = path2CLEO / "libs" / "cleoconstants.hpp"
binpath = path2build / "bin"
sharepath = path2build / "share"
initsupers_filename_1 = sharepath / "shima2009_dimlessSDsinit_1.dat"
initsupers_filename_2 = sharepath / "shima2009_dimlessSDsinit_2.dat"
grid_filename = sharepath / "shima2009_dimlessGBxboundaries.dat"

# path and file names for plotting results
setupfile = binpath / "shima2009_setup.txt"
dataset = binpath / "shima2009_sol.zarr"

# booleans for [making, saving] initialisation figures
isfigures = [True, True]
savefigpath = binpath

### --- settings for 0-D Model gridbox boundaries --- ###
zgrid = np.asarray([0, 100])
xgrid = np.asarray([0, 100])
ygrid = np.asarray([0, 100])

### --- settings for initial superdroplets --- ###
# settings for superdroplet coordinates
nsupers_1 = {0: 4096}
nsupers_2 = {0: 8192}
coord_params = ["false"]

# settings for superdroplet attributes
dryradius = 1e-16  # all SDs have negligible solute [m]
coord3gen = None  # do not generate superdroplet coords
coord1gen = None
coord2gen = None

# radius distirbution from exponential in droplet volume for setup 1
rspan_1 = [0.62e-6, 6.34e-2]  # max and min range of radii to sample [m]
volexpr0_1 = 30.531e-6  # peak of volume exponential distribution [m]
numconc_1 = 2**23  # total no. conc of real droplets [m^-3]
params_1 = {
    "COLLTSTEP": 1,
    "maxnsupers": nsupers_1[0],
    "initsupers_filename": str(initsupers_filename_1),
}

# radius distirbution from exponential in droplet volume for setup 2
rspan_2 = [0.62e-6, 6.34e-2]  # max and min range of radii to sample [m]
volexpr0_2 = 10.117e-6  # peak of volume exponential distribution [m]
numconc_2 = 27 * 2**23  # total no. conc of real droplets [m^-3]
params_2 = {
    "COLLTSTEP": 0.1,
    "maxnsupers": nsupers_2[0],
    "initsupers_filename": str(initsupers_filename_2),
}

# attribute generators
xiprobdist_1 = attrgens_shima2009.SampleXiShima2009()
radiigen_1 = attrgens_shima2009.SampleRadiiShima2009(
    volexpr0_1, rspan_1
)  # radii are sampled from rspan [m]
xiprobdist_2 = attrgens_shima2009.SampleXiShima2009()
radiigen_2 = attrgens_shima2009.SampleRadiiShima2009(
    volexpr0_2, rspan_2
)  # radii are sampled from rspan [m]
samplevol = rgrid.calc_domainvol(zgrid, xgrid, ygrid)
dryradiigen = rgens.MonoAttrGen(dryradius)

### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###

### ---------------------------------------------------------------- ###
### ------------------- BINARY FILES GENERATION--------------------- ###
### ---------------------------------------------------------------- ###
### --- ensure build, share and bin directories exist --- ###
if path2CLEO == path2build:
    raise ValueError("build directory cannot be CLEO")
else:
    path2build.mkdir(exist_ok=True)
    sharepath.mkdir(exist_ok=True)
    binpath.mkdir(exist_ok=True)
    if isfigures[1]:
        savefigpath.mkdir(exist_ok=True)

### --- delete any existing initial conditions --- ###
shutil.rmtree(grid_filename, ignore_errors=True)
shutil.rmtree(initsupers_filename_1, ignore_errors=True)
shutil.rmtree(initsupers_filename_2, ignore_errors=True)

### ----- write gridbox boundaries binary ----- ###
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


### ----- write initial superdroplets binary ----- ###
def initial_conditions_for_setup(
    initsupers_filename, nsupers, radiigen, xiprobdist, numconc, savelabel
):
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
        savelabel=savelabel,
    )


if "golovin" in kernels or "long1" in kernels:
    initial_conditions_for_setup(
        initsupers_filename_1, nsupers_1, radiigen_1, xiprobdist_1, numconc_1, "_1"
    )
if "long2" in kernels:
    initial_conditions_for_setup(
        initsupers_filename_2, nsupers_2, radiigen_2, xiprobdist_2, numconc_2, "_2"
    )
### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###


def run_exectuable(path2build, dataset, executable, config_filename):
    """delete existing dataset, then run exectuable with given config file"""
    os.chdir(path2build)
    subprocess.run(["pwd"])
    shutil.rmtree(dataset, ignore_errors=True)  # delete any existing dataset
    print("Executable: " + str(executable))
    print("Config file: " + str(config_filename))
    subprocess.run([executable, config_filename])


def plot_results(
    grid_filename,
    setupfile,
    dataset,
    numconc,
    volexpr0,
    t2plts,
    smoothsigconst,
    xlims,
    plotwitherr,
    withgol,
    savename,
):
    # read in constants and intial setup from setup .txt file
    config = pysetuptxt.get_config(setupfile, nattrs=3, isprint=True)
    consts = pysetuptxt.get_consts(setupfile, isprint=True)
    gbxs = pygbxsdat.get_gridboxes(grid_filename, consts["COORD0"], isprint=True)

    time = pyzarr.get_time(dataset).secs
    sddata = pyzarr.get_supers(dataset, consts)

    # make and save plot
    savename = savefigpath / savename
    smoothsig = smoothsigconst * (
        config["maxnsupers"] ** (-1 / 5)
    )  # = ~0.2 for guassian smoothing
    fig, ax = shima2009fig.plot_validation_figure(
        plotwitherr,
        time,
        sddata,
        t2plts,
        gbxs["domainvol"],
        numconc,
        volexpr0,
        smoothsig,
        xlims=xlims,
        savename=savename,
        withgol=withgol,
    )


if "golovin" in kernels:
    ### ------------------------------------------------------------ ###
    ### -------------------- RUN CLEO EXECUTABLE ------------------- ###
    ### ------------------------------------------------------------ ###
    editconfigfile.edit_config_params(config_filename, params_1)
    executable = (
        path2build / "examples" / "boxmodelcollisions" / "golovin" / "src" / "golcolls"
    )
    run_exectuable(path2build, dataset, executable, config_filename)
    ### ------------------------------------------------------------ ###
    ### ------------------------------------------------------------ ###

    ### ------------------------------------------------------------ ###
    ### ----------------------- PLOT RESULTS ----------------------- ###
    ### ------------------------------------------------------------ ###
    t2plts = [0, 1200, 2400, 3600]
    smoothsigconst = 0.62
    xlims = [10, 5000]
    plotwitherr = True
    withgol = True
    savename = "golovin_validation.png"
    plot_results(
        grid_filename,
        setupfile,
        dataset,
        numconc_1,
        volexpr0_1,
        t2plts,
        smoothsigconst,
        xlims,
        plotwitherr,
        withgol,
        savename,
    )
    ### ------------------------------------------------------------ ###
    ### ------------------------------------------------------------ ###

if "long1" in kernels:
    ### ------------------------------------------------------------ ###
    ### -------------------- RUN CLEO EXECUTABLE ------------------- ###
    ### ------------------------------------------------------------ ###
    editconfigfile.edit_config_params(config_filename, params_1)
    executable = (
        path2build / "examples" / "boxmodelcollisions" / "long" / "src" / "longcolls"
    )
    run_exectuable(path2build, dataset, executable, config_filename)
    ### ------------------------------------------------------------ ###
    ### ------------------------------------------------------------ ###

    ### ------------------------------------------------------------ ###
    ### ----------------------- PLOT RESULTS ----------------------- ###
    ### ------------------------------------------------------------ ###
    t2plts = [0, 600, 1200, 1800]
    smoothsigconst = 0.62
    xlims = [10, 5000]
    plotwitherr = False
    withgol = False
    savename = "long_validation_1.png"
    plot_results(
        grid_filename,
        setupfile,
        dataset,
        numconc_1,
        volexpr0_1,
        t2plts,
        smoothsigconst,
        xlims,
        plotwitherr,
        withgol,
        savename,
    )
    ### ------------------------------------------------------------ ###
    ### ------------------------------------------------------------ ###

if "long2" in kernels:
    ### ------------------------------------------------------------ ###
    ### -------------------- RUN CLEO EXECUTABLE ------------------- ###
    ### ------------------------------------------------------------ ###
    editconfigfile.edit_config_params(config_filename, params_2)
    executable = (
        path2build / "examples" / "boxmodelcollisions" / "long" / "src" / "longcolls"
    )
    run_exectuable(path2build, dataset, executable, config_filename)
    ### ------------------------------------------------------------ ###
    ### ------------------------------------------------------------ ###

    ### ------------------------------------------------------------ ###
    ### ----------------------- PLOT RESULTS ----------------------- ###
    ### ------------------------------------------------------------ ###
    t2plts = [0, 1200, 1800, 2400, 3600]
    smoothsigconst = 1.0
    xlims = [1, 5000]
    plotwitherr = False
    withgol = False
    savename = "long_validation_2.png"
    plot_results(
        grid_filename,
        setupfile,
        dataset,
        numconc_2,
        volexpr0_2,
        t2plts,
        smoothsigconst,
        xlims,
        plotwitherr,
        withgol,
        savename,
    )
    ### ------------------------------------------------------------ ###
    ### ------------------------------------------------------------ ###
