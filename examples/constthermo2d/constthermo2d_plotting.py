"""
Copyright (c) 2025 MPI-M, Clara Bayley


----- CLEO -----
File: constthermo2d_plotting.py
Project: constthermo2d
Created Date: Friday 22nd August 2025
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
Script plots precipitation example given 2-D flow field and constant
thermodynamics read from a file.
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
        default="/home/m/m300950/CLEO/build_const2d/bin",
    )
    parser.add_argument(
        "--grid_filename",
        type=Path,
        help="Absolute path to gridbox boundaries file",
        default="/home/m/m300950/CLEO/build_const2d/share/const2d_dimlessGBxboundaries.dat",
    )
    parser.add_argument(
        "--setupfile",
        type=Path,
        help="Absolute path to setup file",
        default="/home/m/m300950/CLEO/build_const2d/bin/const2d_setup.txt",
    )
    parser.add_argument(
        "--dataset",
        type=Path,
        help="Absolute path to dataset",
        default="/home/m/m300950/CLEO/build_const2d/bin/const2d_sol.zarr",
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
    from matplotlib.colors import LogNorm, Normalize

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
    savename = savefigpath / "const2d_totnsupers.png"
    pltmoms.plot_totnsupers(time, totnsupers, savename=savename)
    plt.show()

    savename = savefigpath / "const2d_domainmassmoms.png"
    pltmoms.plot_domainmassmoments(time, massmoms, savename=savename)
    plt.show()

    nsample = 500
    savename = savefigpath / "const2d_randomsample.png"
    pltsds.plot_randomsample_superdrops(time, superdrops, nsample, savename=savename)
    plt.show()

    savename = savefigpath / "const2d_motion2d.png"
    superdrops.attach_time(time.secs, "s", do_reshape=True, var4reshape="sdId")
    pltsds.plot_randomsample_superdrops_2dmotion(
        superdrops,
        nsample,
        savename=savename,
        arrows=False,
        cmap_var=["viridis", "time"],
    )
    plt.show()

    ### ----- plot 1-D .gif animations ----- ###
    def horizontal_average(data4d):
        """avg 4-D data with dims [time, y, x, z]
        over x and y dimensions"""
        return np.mean(data4d, axis=(1, 2))

    nframes = len(time.mins)
    norm = np.sum(gbxs["gbxvols"], axis=0)[None, None, :, :] * 1e6  # volume [cm^3]
    mom2ani = horizontal_average(massmoms.mom0 / norm)
    xlims = [0, np.amax(mom2ani)]
    xlabel = "mean number concentration /cm$^{-3}$"
    savename = savefigpath / "const2d_numconc1d"
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

    ### ----- plot 2-D .gif animations ----- ###
    nframes = len(time.mins)
    mom2ani = np.sum(massmoms.nsupers, axis=1)  # sum over y dimension
    cmap = "plasma_r"
    cmapnorm = Normalize(vmin=1, vmax=20)
    cbarlabel = "number of superdroplets per gridbox"
    savename = savefigpath / "const2d_nsupers2d"
    animations.animate2dcmap(
        gbxs,
        mom2ani,
        time.mins,
        nframes,
        cbarlabel=cbarlabel,
        cmapnorm=cmapnorm,
        cmap=cmap,
        saveani=True,
        savename=savename,
        fps=10,
    )

    nframes = len(time.mins)
    mom2ani = np.sum(massmoms.mom1, axis=1)  # sum over y dimension
    norm = np.sum(gbxs["gbxvols"], axis=0)[
        None, :, :
    ]  # sum over y dimension and add time dimension for broadcasting [m^3]
    mom2ani = mom2ani / norm
    cmap = "bone_r"
    cmapnorm = LogNorm(vmin=1e-6, vmax=1e2)
    cbarlabel = "mass concentration /g m$^{-3}$"
    savename = savefigpath / "const2d_massconc2d"
    animations.animate2dcmap(
        gbxs,
        mom2ani,
        time.mins,
        nframes,
        cbarlabel=cbarlabel,
        cmapnorm=cmapnorm,
        cmap=cmap,
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
