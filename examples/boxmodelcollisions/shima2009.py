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
import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

path2CLEO = sys.argv[1]
path2build = sys.argv[2]
configfile = sys.argv[3]
kernels = sys.argv[4:]

sys.path.append(path2CLEO)  # for imports from pySD package
sys.path.append(
    path2CLEO + "/examples/exampleplotting/"
)  # for imports from example plotting package

import attrgens_shima2009
from plotssrc import shima2009fig
from pySD import editconfigfile
from pySD.sdmout_src import pyzarr, pysetuptxt, pygbxsdat
from pySD.initsuperdropsbinary_src import rgens, probdists, attrsgen
from pySD.initsuperdropsbinary_src import create_initsuperdrops as csupers
from pySD.initsuperdropsbinary_src import read_initsuperdrops as rsupers
from pySD.gbxboundariesbinary_src import read_gbxboundaries as rgrid
from pySD.gbxboundariesbinary_src import create_gbxboundaries as cgrid

### ---------------------------------------------------------------- ###
### ----------------------- INPUT PARAMETERS ----------------------- ###
### ---------------------------------------------------------------- ###
### --- essential paths and filenames --- ###
# path and filenames for creating initial SD conditions
constsfile = path2CLEO + "/libs/cleoconstants.hpp"
binpath = path2build + "/bin/"
sharepath = path2build + "/share/"
initSDsfile_1 = sharepath + "shima2009_dimlessSDsinit_1.dat"
initSDsfile_2 = sharepath + "shima2009_dimlessSDsinit_2.dat"
gridfile = sharepath + "shima2009_dimlessGBxboundaries.dat"

# path and file names for plotting results
setupfile = binpath + "shima2009_setup.txt"
dataset = binpath + "shima2009_sol.zarr"

# booleans for [making, saving] initialisation figures
isfigures = [True, True]
savefigpath = path2build + "/bin/"  # directory for saving figures

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
    "initsupers_filename": initSDsfile_1,
}

# radius distirbution from exponential in droplet volume for setup 2
rspan_2 = [5e-7, 5e-5]  # max and min range of radii to sample [m]
volexpr0_2 = 20.117e-6  # peak of volume exponential distribution [m]
numconc_2 = (3 / 2) ** 3 * 2**23  # total no. conc of real droplets [m^-3]
params_2 = {
    "COLLTSTEP": 0.1,
    "maxnsupers": nsupers_2[0],
    "initsupers_filename": initSDsfile_2,
}

# attribute generators
xiprobdist_1 = attrgens_shima2009.SampleXiShima2009()
radiigen_1 = attrgens_shima2009.SampleRadiiShima2009(
    volexpr0_1, rspan_1
)  # radii are sampled from rspan [m]
xiprobdist_2 = probdists.VolExponential(volexpr0_2, rspan_2)
radiigen_2 = rgens.SampleLog10RadiiGen(rspan_2)
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
    Path(path2build).mkdir(exist_ok=True)
    Path(sharepath).mkdir(exist_ok=True)
    Path(binpath).mkdir(exist_ok=True)
    if isfigures[1]:
        Path(savefigpath).mkdir(exist_ok=True)

### --- delete any existing initial conditions --- ###
os.system("rm " + gridfile)
os.system("rm " + initSDsfile_1)
os.system("rm " + initSDsfile_2)

### ----- write gridbox boundaries binary ----- ###
cgrid.write_gridboxboundaries_binary(gridfile, zgrid, xgrid, ygrid, constsfile)
rgrid.print_domain_info(constsfile, gridfile)
### show (and save) plots of binary file data
if isfigures[0]:
    rgrid.plot_gridboxboundaries(constsfile, gridfile, savefigpath, isfigures[1])
    plt.close()


### ----- write initial superdroplets binary ----- ###
def initial_conditions_for_setup(
    initSDsfile, nsupers, radiigen, xiprobdist, numconc, savelabel
):
    initattrsgen = attrsgen.AttrsGenerator(
        radiigen, dryradiigen, xiprobdist, coord3gen, coord1gen, coord2gen
    )
    csupers.write_initsuperdrops_binary(
        initSDsfile, initattrsgen, configfile, constsfile, gridfile, nsupers, numconc
    )
    rsupers.print_initSDs_infos(initSDsfile, configfile, constsfile, gridfile)

    ### show (and save) plots of binary file data
    if isfigures[0]:
        rsupers.plot_initGBxs_distribs(
            configfile,
            constsfile,
            initSDsfile,
            gridfile,
            savefigpath,
            isfigures[1],
            "all",
            savelabel=savelabel,
        )
        plt.close()


if "golovin" in kernels or "long1" in kernels:
    initial_conditions_for_setup(
        initSDsfile_1, nsupers_1, radiigen_1, xiprobdist_1, numconc_1, "_1"
    )
if "long2" in kernels:
    initial_conditions_for_setup(
        initSDsfile_2, nsupers_2, radiigen_2, xiprobdist_2, numconc_2, "_2"
    )
### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###


def run_exectuable(path2build, dataset, executable, configfile):
    """delete existing dataset, the run exectuable with given config file"""
    os.chdir(path2build)
    os.system("pwd")
    os.system("rm -rf " + dataset)  # delete any existing dataset
    print("Executable: " + executable)
    print("Config file: " + configfile)
    os.system(executable + " " + configfile)


def plot_results(
    gridfile,
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
    gbxs = pygbxsdat.get_gridboxes(gridfile, consts["COORD0"], isprint=True)

    time = pyzarr.get_time(dataset).secs
    sddata = pyzarr.get_supers(dataset, consts)

    # make and save plot
    savename = savefigpath + savename
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
    editconfigfile.edit_config_params(configfile, params_1)
    executable = path2build + "/examples/boxmodelcollisions/golovin/src/golcolls"
    run_exectuable(path2build, dataset, executable, configfile)
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
        gridfile,
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
    editconfigfile.edit_config_params(configfile, params_1)
    executable = path2build + "/examples/boxmodelcollisions/long/src/longcolls"
    run_exectuable(path2build, dataset, executable, configfile)
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
        gridfile,
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
    editconfigfile.edit_config_params(configfile, params_2)
    executable = path2build + "/examples/boxmodelcollisions/long/src/longcolls"
    run_exectuable(path2build, dataset, executable, configfile)
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
        gridfile,
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
