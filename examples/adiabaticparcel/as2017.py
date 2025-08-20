"""
Copyright (c) 2024 MPI-M, Clara Bayley


----- CLEO -----
File: as2017.py
Project: adiabaticparcel
Created Date: Friday 17th November 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
Script generate input files, runs CLEO adia0d executable to create data and
then creates plots for adiabatic parcel example similar to Figure 5 of "On
the CCN (de)activation nonlinearities" S. Arabas and S. Shima 2017 to show
example of adaibatic parcel expansion and contraction.
Note: SD(M) = superdroplet (model)
"""

import os
import shutil
import subprocess
import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

path2CLEO = Path(sys.argv[1])
path2build = Path(sys.argv[2])
config_filename = Path(sys.argv[3])

sys.path.append(str(path2CLEO))  # imports from pySD
sys.path.append(
    str(path2CLEO / "examples" / "exampleplotting")
)  # imports from example plots package

from plotssrc import as2017fig
from pySD import editconfigfile, geninitconds
from pySD.initsuperdropsbinary_src import rgens, dryrgens, probdists, attrsgen
from pySD.sdmout_src import pyzarr, pysetuptxt, pygbxsdat

############### INPUTS ##################
# path and filenames for creating SD initial conditions and for running model
constants_filename = path2CLEO / "libs" / "cleoconstants.hpp"
binpath = path2build / "bin"
sharepath = path2build / "share"
initsupers_filename = sharepath / "as2017_dimlessSDsinit.dat"
grid_filename = sharepath / "as2017_dimlessGBxboundaries.dat"

# booleans for [showing, saving] initialisation figures
isfigures = [False, True]
savefigpath = binpath

# settings for 0D Model (number of SD and grid coordinates)
nsupers = {0: 64}
coord_params = ["false"]
zgrid = np.asarray([0, 100])
xgrid = np.asarray([0, 100])
ygrid = np.asarray([0, 100])

# settings for monodisperse droplet radii
# [m^-3] total no. concentration of droplets
numconcs = [500e6, 500e6, 50e6]
monors = [0.05e-6, 0.1e-6, 0.1e-6]

# do not generate superdroplet coords
coord3gen = None
coord1gen = None
coord2gen = None

# parameters to edit in model configuration and plotting
params1 = {
    "W_avg": 1,
    "TAU_half": 150,
    "T_END": 300,
    "COUPLTSTEP": 1,
    "OBSTSTEP": 2,
}
params2 = {
    "W_avg": 0.5,
    "TAU_half": 300,
    "T_END": 600,
    "COUPLTSTEP": 1,
    "OBSTSTEP": 2,
}
params3 = {
    "W_avg": 0.002,
    "TAU_half": 75000,
    "T_END": 150000,
    "COUPLTSTEP": 3,
    "OBSTSTEP": 750,
}
paramslist = [params1, params2, params3]
lwdths = [2, 1, 0.5]


def displacement(time, w_avg, thalf):
    """displacement z given velocity, w, is sinusoidal
    profile: w = w_avg * pi/2 * np.sin(np.pi * t/thalf)
    where wmax = pi/2*w_avg and tauhalf = thalf/pi."""

    zmax = w_avg / 2 * thalf
    z = zmax * (1 - np.cos(np.pi * time / thalf))
    return z


############### RUN EXAMPLE ##################
###  delete any existing datasets
for run_num in range(len(monors) * len(paramslist)):
    dataset = "as2017_sol" + str(run_num) + ".zarr"
    shutil.rmtree(binpath / dataset, ignore_errors=True)

### ensure build, share and bin directories exist
if path2CLEO == path2build:
    raise ValueError("build directory cannot be CLEO")
else:
    path2build.mkdir(exist_ok=True)
    sharepath.mkdir(exist_ok=True)
    binpath.mkdir(exist_ok=True)
    if isfigures[1]:
        savefigpath.mkdir(exist_ok=True)

### create file (and plot) for gridbox boundaries
shutil.rmtree(grid_filename, ignore_errors=True)
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

runnum = 0
for i in range(len(monors)):
    ### create file (and plots) for initial SDs conditions
    monor, numconc = monors[i], numconcs[i]
    # all SDs have the same dryradius = monor [m]
    radiigen = rgens.MonoAttrGen(monor)
    dryradiigen = dryrgens.ScaledRadiiGen(1.0)
    # monodisperse droplet radii probability distribution
    xiprobdist = probdists.DiracDelta(monor)

    initattrsgen = attrsgen.AttrsGenerator(
        radiigen, dryradiigen, xiprobdist, coord3gen, coord1gen, coord2gen
    )
    shutil.rmtree(initsupers_filename, ignore_errors=True)
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

    fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(5, 16))
    for params, lwdth in zip(paramslist, lwdths):
        ### edit relevant setup file parameters
        zarrbasedir = "as2017_sol" + str(runnum) + ".zarr"
        params["zarrbasedir"] = str(binpath / zarrbasedir)
        params["setup_filename"] = str(binpath / "as2017_setup.txt")
        editconfigfile.edit_config_params(config_filename, params)

        ### delete any existing dataset
        shutil.rmtree(params["zarrbasedir"], ignore_errors=True)
        shutil.rmtree(params["setup_filename"], ignore_errors=True)

        ### run model
        os.chdir(path2build)
        executable = path2build / "examples" / "adiabaticparcel" / "src" / "adia0d"
        print("Executable: " + str(executable))
        print("Config file: " + str(config_filename))
        subprocess.run([executable, config_filename])

        ### load results
        setupfile = binpath / "as2017_setup.txt"
        dataset = binpath / Path("as2017_sol" + str(runnum) + ".zarr")

        ### read in constants and intial setup from setup .txt file
        config = pysetuptxt.get_config(setupfile, nattrs=3, isprint=True)
        consts = pysetuptxt.get_consts(setupfile, isprint=True)
        gbxs = pygbxsdat.get_gridboxes(grid_filename, consts["COORD0"], isprint=True)

        # read in output Xarray data
        thermo = pyzarr.get_thermodata(dataset, config["ntime"], gbxs["ndims"], consts)
        supersat = thermo.supersaturation()
        time = pyzarr.get_time(dataset).secs
        superdrops = pyzarr.get_supers(dataset, consts)
        zprof = displacement(time, config["W_avg"], config["TAU_half"])

        attrs = ["radius", "xi", "msol"]
        sd0 = superdrops.sample("sdId", sample_values=0, variables2sample=attrs)
        radius = sd0["radius"][0]
        msol_initial = sd0["msol"][0][0]
        numconc = np.sum(superdrops["xi"][0]) / gbxs["domainvol"] / 1e6  # [/cm^3]

        ### plot results
        wlab = "<w> = {:.1f}".format(config["W_avg"] * 100) + "cm s$^{-1}$"
        as2017fig.condensation_validation_subplots(
            axs, time, radius, supersat[:, 0, 0, 0], zprof, lwdth=lwdth, lab=wlab
        )

        runnum += 1

    ### save figure
    as2017fig.plot_kohlercurve_with_criticalpoints(
        axs[1],
        radius,
        msol_initial,
        thermo.temp[0, 0, 0, 0],
        superdrops.IONIC(),
        superdrops.MR_SOL(),
    )

    textlab = (
        "N = "
        + str(numconc)
        + "cm$^{-3}$\n"
        + "r$_{dry}$ = "
        + "{:.2g}\u03BCm\n".format(radius[0])
    )
    axs[0].legend(loc="lower right", fontsize=10)
    axs[1].legend(loc="upper left")
    axs[0].text(0.03, 0.85, textlab, transform=axs[0].transAxes)

    axs[0].set_xlim([-1, 1])
    for ax in axs[1:]:
        ax.set_xlim([0.125, 10])
        ax.set_xscale("log")
    axs[0].set_ylim([0, 150])
    axs[1].set_ylim([-1, 1])
    axs[2].set_ylim([5, 75])

    fig.tight_layout()

    savename = savefigpath / Path("as2017fig_" + str(i) + ".png")
    fig.savefig(savename, dpi=400, bbox_inches="tight", facecolor="w", format="png")
    print("Figure .png saved as: " + str(savename))
    plt.show()
