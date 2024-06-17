"""
Copyright (c) 2024 MPI-M, Clara Bayley


----- CLEO -----
File: breakup.py
Project: boxmodelcollisions
Created Date: Friday 14th June 2024
Author: Clara Bayley (CB)
Additional Contributors:
-----
Last Modified: Monday 17th June 2024
Modified By: CB
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
import sys
import numpy as np
import awkward as ak
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

from plotssrc import shima2009fig
from pySD import editconfigfile
from pySD.sdmout_src import pyzarr, pysetuptxt, pygbxsdat, sdtracing
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
initSDsfile = sharepath + "breakup_dimlessSDsinit.dat"
gridfile = sharepath + "breakup_dimlessGBxboundaries.dat"

# booleans for [making, saving] initialisation figures
isfigures = [True, True]
savefigpath = path2build + "/bin/"  # directory for saving figures

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
    Path(path2build).mkdir(exist_ok=True)
    Path(sharepath).mkdir(exist_ok=True)
    Path(binpath).mkdir(exist_ok=True)
    if isfigures[1]:
        Path(savefigpath).mkdir(exist_ok=True)

### --- delete any existing initial conditions --- ###
os.system("rm " + gridfile)
os.system("rm " + initSDsfile)

### ----- write gridbox boundaries binary ----- ###
cgrid.write_gridboxboundaries_binary(gridfile, zgrid, xgrid, ygrid, constsfile)
rgrid.print_domain_info(constsfile, gridfile)

### ----- write initial superdroplets binary ----- ###
initattrsgen = attrsgen.AttrsGenerator(
    radiigen, dryradiigen, xiprobdist, coord3gen, coord1gen, coord2gen
)
csupers.write_initsuperdrops_binary(
    initSDsfile, initattrsgen, configfile, constsfile, gridfile, nsupers, numconc
)
rsupers.print_initSDs_infos(initSDsfile, configfile, constsfile, gridfile)

### show (and save) plots of binary file data
if isfigures[0]:
    rgrid.plot_gridboxboundaries(constsfile, gridfile, savefigpath, isfigures[1])
    rsupers.plot_initGBxs_distribs(
        configfile,
        constsfile,
        initSDsfile,
        gridfile,
        savefigpath,
        isfigures[1],
        "all",
    )
    plt.close()
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
        path2build
        + "/examples/boxmodelcollisions/"
        + kernel
        + "/src/"
        + executables[kernel]
    )

    return executable


def get_params(path2build, kernel):
    params = {
        "setup_filename": path2build + "bin/" + kernel + "_setup.txt",
        "stats_filename": path2build + "bin/" + kernel + "_stats.txt",
        "zarrbasedir": path2build + "bin/" + kernel + "_sol.zarr",
    }

    return params


def run_exectuable(path2build, kernel, configfile):
    """delete existing dataset, the run exectuable with given config file"""
    params = get_params(path2build, kernel)
    editconfigfile.edit_config_params(configfile, params)

    executable = get_executable(path2build, kernel)
    os.chdir(path2build)
    os.system("pwd")
    os.system("rm -rf " + params["zarrbasedir"])  # delete any existing dataset
    print("Executable: " + executable)
    print("Config file: " + configfile)
    os.system(executable + " " + configfile)


def get_kernel_results(path2build, kernel):
    """read in time, sddata and setup for given kernel"""
    params = get_params(path2build, kernel)
    setupfile = params["setup_filename"]
    dataset = params["zarrbasedir"]

    # read in constants and intial setup from setup .txt file
    config = pysetuptxt.get_config(setupfile, nattrs=3, isprint=True)
    consts = pysetuptxt.get_consts(setupfile, isprint=True)

    # read in seom data from dataset
    time = pyzarr.get_time(dataset).secs
    sddata = pyzarr.get_supers(dataset, consts)

    return config, consts, time, sddata


def plot_onekernel_results(
    gridfile,
    path2build,
    kernel,
    numconc,
    volexpr0,
    xlims,
    t2plts,
    savename,
):
    # read in data
    config, consts, time, sddata = get_kernel_results(path2build, kernel)
    gbxs = pygbxsdat.get_gridboxes(gridfile, consts["COORD0"], isprint=True)

    # make and save plot
    smoothsigconst = 0.62
    smoothsig = smoothsigconst * (
        config["maxnsupers"] ** (-1 / 5)
    )  # = ~0.2 for guassian smoothing
    plotwitherr = False
    withgol = False
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


def plot_allkernels_results(gridfile, path2build, kernels, xlims, t2plts, savename):
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
        config, consts, time, sddata = get_kernel_results(path2build, kernel)
        domainvol = pygbxsdat.get_gridboxes(gridfile, consts["COORD0"], isprint=True)[
            "domainvol"
        ]
        smoothsig = 0.62 * (
            config["maxnsupers"] ** (-1 / 5)
        )  # = ~0.2 for guassian smoothing

        attrs2sel = ["radius", "xi"]
        selsddata = sdtracing.attributes_at_times(sddata, time, t2plts, attrs2sel)

        nbins = 500
        non_nanradius = ak.nan_to_none(sddata["radius"])
        rspan = [ak.min(non_nanradius), ak.max(non_nanradius)]

        for n in range(len(t2plts)):
            radius = selsddata["radius"][n]
            xi = selsddata["xi"][n]
            hist, hcens = shima2009fig.calc_massdens_distrib(
                rspan, nbins, domainvol, xi, radius, sddata, smoothsig
            )
            if kernel == "long":
                ind = np.argmin(abs(time - t2plts[n]))
                tlab = "t = {:.2f}s".format(time[ind])
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
        print("Figure .png saved as: " + savename)
    plt.show()

    return fig, ax


### ------------------------------------------------------------ ###
### ---------- RUN CLEO EXECUTABLES FOR EACH KERNEL ------------ ###
### ------------------------------------------------------------ ###
for kernel in kernels:
    run_exectuable(path2build, kernel, configfile)
### ------------------------------------------------------------ ###
### ------------------------------------------------------------ ###

### ------------------------------------------------------------ ###
### ----------------------- PLOT RESULTS ----------------------- ###
### ------------------------------------------------------------ ###
for kernel in kernels:
    t2plts = [0, 600, 1200, 1800, 2400]
    xlims = [10, 5000]
    savename = savefigpath + kernel + "_validation.png"
    plot_onekernel_results(
        gridfile,
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
savename = savefigpath + "breakup_validation.png"
plot_allkernels_results(gridfile, path2build, kernels, xlims, t2plts, savename)
### ------------------------------------------------------------ ###
### ------------------------------------------------------------ ###
