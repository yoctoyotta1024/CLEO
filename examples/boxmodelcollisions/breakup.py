"""
Copyright (c) 2024 MPI-M, Clara Bayley


----- CLEO -----
File: breakup.py
Project: boxmodelcollisions
Created Date: Friday 14th June 2024
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
Script generates input files, runs CLEO 0-D box model executables for
collisions with selected collision kernels with breakup (e.g. Low and Lists's)
to create data. Then plots results analgous to Shima et al. 2009 Fig. 2(b)
"""

import os
import shutil
import subprocess
import sys
import numpy as np
import awkward as ak
import matplotlib.pyplot as plt
from pathlib import Path

path2CLEO = Path(sys.argv[1])
path2build = Path(sys.argv[2])
config_filename = Path(sys.argv[3])
kernels = sys.argv[4:]

sys.path.append(str(path2CLEO))  # imports from pySD
sys.path.append(
    str(path2CLEO / "examples" / "exampleplotting")
)  # imports from example plots package


from plotssrc import shima2009fig
from pySD import editconfigfile, geninitconds
from pySD.sdmout_src import pyzarr, pysetuptxt, pygbxsdat
from pySD.initsuperdropsbinary_src import rgens, probdists, attrsgen
from pySD.gbxboundariesbinary_src import read_gbxboundaries as rgrid

### ---------------------------------------------------------------- ###
### ----------------------- INPUT PARAMETERS ----------------------- ###
### ---------------------------------------------------------------- ###
### --- essential paths and filenames --- ###
# path and filenames for creating initial SD conditions
constants_filename = path2CLEO / "libs" / "cleoconstants.hpp"
binpath = path2build / "bin"
sharepath = path2build / "share"
initsupers_filename = sharepath / "breakup_dimlessSDsinit.dat"
grid_filename = sharepath / "breakup_dimlessGBxboundaries.dat"

# booleans for [showing, saving] initialisation figures
isfigures = [False, True]
savefigpath = binpath

### --- settings for 0-D Model gridbox boundaries --- ###
zgrid = np.asarray([0, 100])
xgrid = np.asarray([0, 100])
ygrid = np.asarray([0, 100])

### --- settings for initial superdroplets --- ###
# settings for superdroplet coordinates
nsupers = {0: 8192}
coord_params = ["false"]

# settings for superdroplet attributes
dryradius = 1e-16  # all SDs have negligible solute [m]
coord3gen = None  # do not generate superdroplet coords
coord1gen = None
coord2gen = None

# radius distirbution from exponential in droplet volume for setup 1
rspan = [1e-7, 9e-5]  # max and min range of radii to sample [m]
volexpr0 = 30.531e-6  # peak of volume exponential distribution [m]
numconc = 2 ** (23)  # total no. conc of real droplets [m^-3]

# attribute generators
xiprobdist = probdists.VolExponential(volexpr0, rspan)
radiigen = rgens.SampleLog10RadiiGen(rspan)  # radii are sampled from rspan [m]
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
shutil.rmtree(initsupers_filename, ignore_errors=True)

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
### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###


def get_executable(path2build, kernel):
    executables = {
        "long": "longcolls",
        "lowlist": "lowlistcolls",
        "szakallurbich": "szakallurbichcolls",
        "testikstraub": "testikstraubcolls",
    }
    executable = (
        path2build / "examples" / "boxmodelcollisions" / "src" / executables[kernel]
    )

    return executable


def get_params(path2build, kernel):
    params = {
        "setup_filename": str(path2build / "bin" / Path(kernel + "_setup.txt")),
        "zarrbasedir": str(path2build / "bin" / Path(kernel + "_sol.zarr")),
    }

    return params


def run_exectuable(path2build, kernel, config_filename):
    """delete existing dataset, the run exectuable with given config file"""
    params = get_params(path2build, kernel)
    editconfigfile.edit_config_params(config_filename, params)

    executable = get_executable(path2build, kernel)
    os.chdir(path2build)
    subprocess.run(["pwd"])
    shutil.rmtree(
        params["zarrbasedir"], ignore_errors=True
    )  # delete any existing dataset
    print("Executable: " + str(executable))
    print("Config file: " + str(config_filename))
    subprocess.run([executable, config_filename])


def get_kernel_results(path2build, kernel):
    """read in time, superdrops and setup for given kernel"""
    params = get_params(path2build, kernel)
    setupfile = params["setup_filename"]
    dataset = params["zarrbasedir"]

    # read in constants and intial setup from setup .txt file
    config = pysetuptxt.get_config(setupfile, nattrs=3, isprint=True)
    consts = pysetuptxt.get_consts(setupfile, isprint=True)

    # read in seom data from dataset
    time = pyzarr.get_time(dataset)
    superdrops = pyzarr.get_supers(dataset, consts)

    return config, consts, time, superdrops


def plot_onekernel_results(
    grid_filename,
    path2build,
    kernel,
    numconc,
    volexpr0,
    xlims,
    t2plts,
    savename,
):
    # read in data
    config, consts, time, superdrops = get_kernel_results(path2build, kernel)
    gbxs = pygbxsdat.get_gridboxes(grid_filename, consts["COORD0"], isprint=True)

    # make and save plot
    smoothsigconst = 0.62
    smoothsig = smoothsigconst * (
        config["maxnsupers"] ** (-1 / 5)
    )  # = ~0.2 for guassian smoothing
    plotwitherr = False
    withgol = False
    shima2009fig.plot_validation_figure(
        plotwitherr,
        time,
        superdrops,
        t2plts,
        gbxs["domainvol"],
        numconc,
        volexpr0,
        smoothsig,
        xlims=xlims,
        savename=savename,
        withgol=withgol,
    )
    plt.show()


def plot_allkernels_results(
    grid_filename, path2build, kernels, xlims, t2plts, savename
):
    def blank_axis(ax, xlims, ylims):
        ax2.set_xlim(xlims)
        ax2.set_ylim(ylims)
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_xticks([])  # Remove x-axis tick marks
        ax.set_yticks([])  # Remove y-axis tick marks
        ax.set_xticklabels([])  # Remove x-axis tick labels
        ax.set_yticklabels([])  # Remove y-axis tick labels

    styles = {
        "long": ["Long (coal only)", "grey", 0.8],
        "lowlist": ["Low & List", "blue", 1.0],
        "szakallurbich": ["Szakall & Urbich", "green", 1.0],
        "testikstraub": ["Testik + Straub", "purple", 1.0],
    }
    linestyles = ["dotted", (0, (3, 1, 1, 1)), "dashed", "dashdot", "solid"]
    witherr = False

    fig, ax = shima2009fig.setup_validation_figure(witherr, xlims)[0:2]
    ax2 = ax.twinx()
    for kernel in kernels:
        # read in data
        config, consts, time, superdrops = get_kernel_results(path2build, kernel)
        domainvol = pygbxsdat.get_gridboxes(
            grid_filename, consts["COORD0"], isprint=True
        )["domainvol"]
        smoothsig = 0.62 * (
            config["maxnsupers"] ** (-1 / 5)
        )  # = ~0.2 for guassian smoothing

        variables2slice = ["time", "radius", "xi"]
        superdrops_timeslice = superdrops.time_slice(
            t2plts, variables2slice, attach_time=True, time=time.secs, time_units="s"
        )

        nbins = 500
        non_nanradius = ak.nan_to_none(superdrops["radius"])
        rspan = [ak.min(non_nanradius), ak.max(non_nanradius)]

        for n in range(len(t2plts)):
            radius = superdrops_timeslice["radius"][n]
            xi = superdrops_timeslice["xi"][n]
            hist, hcens = shima2009fig.calc_massdens_distrib(
                rspan, nbins, domainvol, xi, radius, superdrops, smoothsig
            )
            if kernel == "long":
                t = superdrops_timeslice.time()[n]
                tlab = "t = {:.2f}".format(t) + superdrops_timeslice.time_units()
                ax2.plot(
                    hcens,
                    hist,
                    label=tlab,
                    color="k",
                    linestyle=linestyles[n],
                    linewidth=0.5,
                )
            if linestyles[n] == "solid":
                ax.plot(
                    hcens,
                    hist,
                    label=styles[kernel][0],
                    color=styles[kernel][1],
                    linestyle=linestyles[n],
                    linewidth=styles[kernel][2],
                )
            else:
                ax.plot(
                    hcens,
                    hist,
                    color=styles[kernel][1],
                    linestyle=linestyles[n],
                    linewidth=styles[kernel][2],
                )

    ax.legend(loc="lower left")
    ax2.legend(loc="lower right")
    blank_axis(ax2, ax.get_xlim(), ax.get_ylim())
    fig.tight_layout()

    if savename != "":
        fig.savefig(savename, dpi=400, bbox_inches="tight", facecolor="w", format="png")
        print("Figure .png saved as: " + str(savename))
    plt.show()

    return fig, ax


### ------------------------------------------------------------ ###
### ---------- RUN CLEO EXECUTABLES FOR EACH KERNEL ------------ ###
### ------------------------------------------------------------ ###
for kernel in kernels:
    run_exectuable(path2build, kernel, config_filename)
### ------------------------------------------------------------ ###
### ------------------------------------------------------------ ###

### ------------------------------------------------------------ ###
### ----------------------- PLOT RESULTS ----------------------- ###
### ------------------------------------------------------------ ###
for kernel in kernels:
    t2plts = [0, 600, 1200, 1800, 2400]
    xlims = [10, 5000]
    savename = savefigpath / Path(kernel + "_validation.png")
    plot_onekernel_results(
        grid_filename,
        path2build,
        kernel,
        numconc,
        volexpr0,
        xlims,
        t2plts,
        savename,
    )


t2plts = [0, 600, 1200, 1800, 2400]
xlims = [10, 5000]
savename = savefigpath / "breakup_validation.png"
plot_allkernels_results(grid_filename, path2build, kernels, xlims, t2plts, savename)
### ------------------------------------------------------------ ###
### ------------------------------------------------------------ ###
