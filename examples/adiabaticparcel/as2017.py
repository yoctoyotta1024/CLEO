"""
Copyright (c) 2024 MPI-M, Clara Bayley


----- CLEO -----
File: as2017.py
Project: adiabaticparcel
Created Date: Friday 17th November 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
Last Modified: Tuesday 7th May 2024
Modified By: CB
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
Script generate input files, runs CLEO adia0D executable to create data and
then creates plots for adiabatic parcel example similar to Figure 5 of "On
the CCN (de)activation nonlinearities" S. Arabas and S. Shima 2017 to show
example of adaibatic parcel expansion and contraction.
Note: SD(M) = superdroplet (model)
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

path2CLEO = sys.argv[1]
path2build = sys.argv[2]
configfile = sys.argv[3]

sys.path.append(path2CLEO)  # for imports from pySD package
sys.path.append(
    path2CLEO + "/examples/exampleplotting/"
)  # for imports from example plotting package

from plotssrc import as2017fig
from pySD import editconfigfile
from pySD.sdmout_src import pyzarr, pysetuptxt, pygbxsdat, sdtracing
from pySD.initsuperdropsbinary_src import rgens, dryrgens, probdists, attrsgen
from pySD.initsuperdropsbinary_src import create_initsuperdrops as csupers
from pySD.initsuperdropsbinary_src import read_initsuperdrops as rsupers
from pySD.gbxboundariesbinary_src import read_gbxboundaries as rgrid
from pySD.gbxboundariesbinary_src import create_gbxboundaries as cgrid

############### INPUTS ##################
# path and filenames for creating SD initial conditions and for running model
constsfile = path2CLEO + "/libs/cleoconstants.hpp"
binpath = path2build + "/bin/"
sharepath = path2build + "/share/"
initSDsfile = sharepath + "as2017_dimlessSDsinit.dat"
gridfile = sharepath + "as2017_dimlessGBxboundaries.dat"

# booleans for [making, saving] initialisation figures
isfigures = [True, True]

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

# volume SD sample occupies (entire domain) [m^3]
samplevol = rgrid.calc_domainvol(zgrid, xgrid, ygrid)
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
    dataset = binpath + "as2017_sol" + str(run_num) + ".zarr"
    os.system("rm -rf " + dataset)

### ensure build, share and bin directories exist
if path2CLEO == path2build:
    raise ValueError("build directory cannot be CLEO")
else:
    Path(path2build).mkdir(exist_ok=True)
    Path(sharepath).mkdir(exist_ok=True)
    Path(binpath).mkdir(exist_ok=True)

### create file (and plot) for gridbox boundaries
os.system("rm " + gridfile)
cgrid.write_gridboxboundaries_binary(gridfile, zgrid, xgrid, ygrid, constsfile)
rgrid.print_domain_info(constsfile, gridfile)
if isfigures[0]:
    rgrid.plot_gridboxboundaries(constsfile, gridfile, binpath, isfigures[1])
plt.close()

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
    os.system("rm " + initSDsfile)
    csupers.write_initsuperdrops_binary(
        initSDsfile, initattrsgen, configfile, constsfile, gridfile, nsupers, numconc
    )
    rsupers.print_initSDs_infos(initSDsfile, configfile, constsfile, gridfile)

    if isfigures[0]:
        rsupers.plot_initGBxs_distribs(
            configfile, constsfile, initSDsfile, gridfile, binpath, isfigures[1], "all"
        )
        plt.close()

    fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(5, 16))
    for params, lwdth in zip(paramslist, lwdths):
        ### edit relevant setup file parameters
        params["zarrbasedir"] = binpath + "as2017_sol" + str(runnum) + ".zarr"
        params["setup_filename"] = binpath + "as2017_setup.txt"
        editconfigfile.edit_config_params(configfile, params)

        ### delete any existing dataset
        os.system("rm -rf " + params["zarrbasedir"])
        os.system("rm " + params["setup_filename"])

        ### run model
        os.chdir(path2build)
        executable = path2build + "/examples/adiabaticparcel/src/adia0D"
        print("Executable: " + executable)
        print("Config file: " + configfile)
        os.system(executable + " " + configfile)

        ### load results
        setupfile = binpath + "as2017_setup.txt"
        dataset = binpath + "as2017_sol" + str(runnum) + ".zarr"

        ### read in constants and intial setup from setup .txt file
        config = pysetuptxt.get_config(setupfile, nattrs=3, isprint=True)
        consts = pysetuptxt.get_consts(setupfile, isprint=True)
        gbxs = pygbxsdat.get_gridboxes(gridfile, consts["COORD0"], isprint=True)

        # read in output Xarray data
        thermo = pyzarr.get_thermodata(dataset, config["ntime"], gbxs["ndims"], consts)
        supersat = thermo.supersaturation()
        time = pyzarr.get_time(dataset).secs
        sddata = pyzarr.get_supers(dataset, consts)
        zprof = displacement(time, config["W_avg"], config["TAU_half"])

        attrs = ["radius", "xi", "msol"]
        sd0 = sdtracing.attributes_for1superdroplet(sddata, 0, attrs)
        numconc = np.sum(sddata["xi"][0]) / gbxs["domainvol"] / 1e6  # [/cm^3]

        ### plot results
        wlab = "<w> = {:.1f}".format(config["W_avg"] * 100) + "cm s$^{-1}$"
        axs = as2017fig.condensation_validation_subplots(
            axs, time, sd0["radius"], supersat[:, 0, 0, 0], zprof, lwdth=lwdth, lab=wlab
        )

        runnum += 1

    ### save figure
    as2017fig.plot_kohlercurve_with_criticalpoints(
        axs[1],
        sd0["radius"],
        sd0["msol"][0],
        thermo.temp[0, 0, 0, 0],
        sddata.IONIC,
        sddata.MR_SOL,
    )

    textlab = (
        "N = "
        + str(numconc)
        + "cm$^{-3}$\n"
        + "r$_{dry}$ = "
        + "{:.2g}\u03BCm\n".format(sd0["radius"][0])
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

    savename = "as2017fig_" + str(i) + ".png"
    fig.savefig(
        binpath + savename, dpi=400, bbox_inches="tight", facecolor="w", format="png"
    )
    print("Figure .png saved as: " + binpath + savename)
    plt.show()
