"""
Copyright (c) 2025 MPI-M, Clara Bayley


----- CLEO -----
File: cuspbifurc_plotting.py
Project: adiabaticparcel
Created Date: Monday 25th August 2025
Author: Clara Bayley (CB)
Additional Contributors:
-----
Last Modified: Monday 25th August 2025
Modified By: CB
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
Script creates plots for adiabatic parcel example similar to Figure 5 of
"On the CCN (de)activation nonlinearities" S. Arabas and S. Shima 2017 to show
example of cusp birfucation for 0D adiabatic parcel expansion and contraction.
"""


# %%
### ------------------------- FUNCTION DEFINITIONS ------------------------- ###
def parse_arguments():
    import argparse
    from pathlib import Path

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--path2CLEO",
        type=Path,
        help="Absolute path to CLEO",
        default="/home/m/m300950/CLEO",
    )
    parser.add_argument(
        "--savefigpath",
        type=Path,
        help="Absolute path to build",
        default="/home/m/m300950/CLEO/build_adia0d/cuspbifurc/bin",
    )
    parser.add_argument(
        "--grid_filename",
        type=Path,
        help="Absolute path to gridbox boundaries file",
        default="/home/m/m300950/CLEO/build_adia0d/cuspbifurc/share/cuspbifurc_dimlessGBxboundaries.dat",
    )
    parser.add_argument(
        "--setupfile",
        type=Path,
        help="Absolute path to setup file",
        default="/home/m/m300950/CLEO/build_adia0d/cuspbifurc/bin/cuspbifurc_setup.txt",
    )
    parser.add_argument(
        "--dataset",
        type=Path,
        help="Absolute path to dataset",
        default="/home/m/m300950/CLEO/build_adia0d/cuspbifurc/bin/cuspbifurc_sol.zarr",
    )
    return parser.parse_args()


def displacement(time, w_avg, thalf):
    """displacement z given velocity, w, is sinusoidal
    profile: w = w_avg * pi/2 * np.sin(np.pi * t/thalf)
    where wmax = pi/2*w_avg and tauhalf = thalf/pi."""
    import numpy as np

    zmax = w_avg / 2 * thalf
    z = zmax * (1 - np.cos(np.pi * time / thalf))
    return z


# %%
### -------------------------------- MAIN ---------------------------------- ###
def main(
    path2CLEO,
    savefigpath,
    grid_filename,
    setupfile,
    dataset,
):
    import sys
    import numpy as np
    import matplotlib.pyplot as plt

    sys.path.append(
        str(path2CLEO / "examples" / "exampleplotting")
    )  # imports from example plots package
    from plotssrc import pltsds, as2017fig
    from cleopy.sdmout_src import pyzarr, pysetuptxt, pygbxsdat

    ### load results
    # read in constants and intial setup from setup .txt file
    config = pysetuptxt.get_config(setupfile, nattrs=3, isprint=True)
    consts = pysetuptxt.get_consts(setupfile, isprint=True)
    gbxs = pygbxsdat.get_gridboxes(grid_filename, consts["COORD0"], isprint=True)

    # read in output Xarray data
    thermo = pyzarr.get_thermodata(dataset, config["ntime"], gbxs["ndims"], consts)
    supersat = thermo.supersaturation()
    time = pyzarr.get_time(dataset).secs
    superdrops = pyzarr.get_supers(dataset, consts)
    zprof = displacement(time, config["W_avg"], config["TAU_half"])

    ### plot results
    # sample drops to plot from whole range of SD ids
    radii = superdrops.sample(
        "sdId", sample_values="all", variables2sample="radius"
    ).radius()
    savename = savefigpath / "cuspbifurc_SDgrowth.png"
    pltsds.individ_radiusgrowths_figure(time, radii, savename=savename)
    plt.show()

    attrs = ["radius", "xi", "msol"]
    sd0 = superdrops.sample("sdId", sample_values=0, variables2sample=attrs)
    radius = np.array(sd0["radius"])  # convert list of awkward arrays to numpy
    msol = np.array(sd0["msol"])  # convert list of awkward arrays to numpy
    numconc = np.sum(superdrops["xi"][0]) / gbxs["domainvol"] / 1e6  # [/cm^3]

    savename2 = savefigpath / "cuspbifurc_validation.png"
    as2017fig.arabas_shima_2017_fig(
        time,
        zprof,
        radius,
        msol,
        thermo.temp[:, 0, 0, 0],
        supersat[:, 0, 0, 0],
        superdrops.IONIC(),
        superdrops.MR_SOL(),
        config["W_avg"],
        numconc,
        savename2,
    )
    plt.show()


# %%
### --------------------------- RUN PROGRAM -------------------------------- ###
if __name__ == "__main__":
    args = parse_arguments()
    main(
        args.path2CLEO,
        args.savefigpath,
        args.grid_filename,
        args.setupfile,
        args.dataset,
    )
