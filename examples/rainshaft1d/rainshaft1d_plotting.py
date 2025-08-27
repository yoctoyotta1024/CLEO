"""
Copyright (c) 2025 MPI-M, Clara Bayley


----- CLEO -----
File: rainshaft1d_plotting.py
Project: rainshaft1d
Created Date: Friday 22nd August 2025
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
Script plots results of example of CLEO executable "rshaft1d" for 1-D rainshaft
with constant thermodynamics read from a file.
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
        default="/home/m/m300950/CLEO/build_rshaft1d/bin",
    )
    parser.add_argument(
        "--grid_filename",
        type=Path,
        help="Absolute path to gridbox boundaries file",
        default="/home/m/m300950/CLEO/build_rshaft1d/share/rshaft1d_dimlessGBxboundaries.dat",
    )
    parser.add_argument(
        "--setupfile",
        type=Path,
        help="Absolute path to setup file",
        default="/home/m/m300950/CLEO/build_rshaft1d/bin/rshaft1d_setup.txt",
    )
    parser.add_argument(
        "--dataset",
        type=Path,
        help="Absolute path to dataset",
        default="/home/m/m300950/CLEO/build_rshaft1d/bin/rshaft1d_sol.zarr",
    )
    return parser.parse_args()


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

    from plotssrc import pltsds, pltmoms, animations
    from cleopy.sdmout_src import pyzarr, pysetuptxt, pygbxsdat

    # read in constants and intial setup from setup .txt file
    config = pysetuptxt.get_config(setupfile, nattrs=3, isprint=True)
    consts = pysetuptxt.get_consts(setupfile, isprint=True)
    gbxs = pygbxsdat.get_gridboxes(grid_filename, consts["COORD0"], isprint=True)

    time = pyzarr.get_time(dataset)
    superdrops = pyzarr.get_supers(dataset, consts)
    totnsupers = pyzarr.get_totnsupers(dataset)
    massmoms = pyzarr.get_massmoms(dataset, config["ntime"], gbxs["ndims"])

    # plot figures
    savename = savefigpath / "rshaft1d_totnsupers.png"
    pltmoms.plot_totnsupers(time, totnsupers, savename=savename)
    plt.show()

    savename = savefigpath / "rshaft1d_domainmassmoms.png"
    pltmoms.plot_domainmassmoments(time, massmoms, savename=savename)
    plt.show()

    nsample = 500
    savename = savefigpath / "rshaft1d_randomsample.png"
    pltsds.plot_randomsample_superdrops(time, superdrops, nsample, savename=savename)
    plt.show()

    ### ----- plot 1-D .gif animations ----- ###
    nframes = len(time.mins)
    mom2ani = np.sum(massmoms.nsupers, axis=(1, 2))
    xlims = [0, np.amax(mom2ani)]
    xlabel = "number of super-droplets"
    savename = savefigpath / "rshaft1d_nsupers1d"
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
        fps=10,
    )

    nframes = len(time.mins)
    norm = gbxs["gbxvols"] * 1e6  # volume [cm^3]
    mom2ani = np.sum(massmoms.mom0 / norm[None, :], axis=(1, 2))
    xlims = [0, np.amax(mom2ani)]
    xlabel = "number concentration /cm$^{-3}$"
    savename = savefigpath / "rshaft1d_numconc1d"
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
        fps=10,
    )

    nframes = len(time.mins)
    norm = gbxs["gbxvols"]  # volume [m^3]
    mom2ani = np.sum(massmoms.mom1 / norm[None, :], axis=(1, 2))
    xlims = [0, np.amax(mom2ani)]
    xlabel = "mass concentration /g m$^{-3}$"
    savename = savefigpath / "rshaft1d_massconc1d"
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
        fps=10,
    )


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
