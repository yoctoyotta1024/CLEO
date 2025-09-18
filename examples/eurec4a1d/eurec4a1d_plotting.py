"""
Copyright (c) 2025 MPI-M, Clara Bayley


----- CLEO -----
File: eurec4a1d_plotting.py
Project: eurec4a1d
Created Date: Thursday 1st January 1970
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
Script plots results of example of EUREC4A 1-D rainshaft (test)
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
        default="/home/m/m300950/CLEO/build_eurec4a1d/bin",
    )
    parser.add_argument(
        "--grid_filename",
        type=Path,
        help="Absolute path to gridbox boundaries file",
        default="/home/m/m300950/CLEO/build_eurec4a1d/share/eurec4a1d_dimlessGBxboundaries.dat",
    )
    parser.add_argument(
        "--setupfile",
        type=Path,
        help="Absolute path to setup file",
        default="/home/m/m300950/CLEO/build_eurec4a1d/bin/eurec4a1d_setup.txt",
    )
    parser.add_argument(
        "--dataset",
        type=Path,
        help="Absolute path to dataset",
        default="/home/m/m300950/CLEO/build_eurec4a1d/bin/eurec4a1d_sol.zarr",
    )
    return parser.parse_args()


# %%
def plot_droplet_distributions(config, gbxs, time, superdrops, savename=""):
    import awkward as ak
    import matplotlib.pyplot as plt
    import numpy as np
    from plotcleo import pltdist

    fig, axs = plt.subplots(nrows=3, ncols=3, figsize=(16, 12), sharex=True)

    t2plts = [time.secs[0], time.secs[10], time.secs[-1]]

    var_for_range = "coord3"
    coordlims_0 = [config["UPPER_COORD3LIM"], np.amax(gbxs["zhalf"])]
    coordlims_1 = [config["LOWER_COORD3LIM"], config["UPPER_COORD3LIM"]]
    coordlims_2 = [np.amin(gbxs["zhalf"]), config["LOWER_COORD3LIM"]]

    rspan = [ak.min(superdrops["radius"]) * 0.9, ak.max(superdrops["radius"]) * 1.1]
    nbins = 100
    smoothsig = False
    perlogR = False
    ylog_nsupers = False
    ylog_num = False
    ylog_mass = True

    for r, coordlims in enumerate([coordlims_0, coordlims_1, coordlims_2]):
        volume = (coordlims[1] - coordlims[0]) * gbxs["domainarea"]
        pltdist.plot_distribs_in_select_range(
            time,
            pltdist.plot_domainnsupers_distribs,
            var_for_range,
            coordlims,
            superdrops,
            t2plts,
            volume,
            rspan,
            nbins,
            smoothsig=smoothsig,
            perlogR=perlogR,
            ylog=ylog_nsupers,
            fig_ax=[fig, axs[r, 0]],
        )
        pltdist.plot_distribs_in_select_range(
            time,
            pltdist.plot_domainnumconc_distribs,
            var_for_range,
            coordlims,
            superdrops,
            t2plts,
            volume,
            rspan,
            nbins,
            smoothsig=smoothsig,
            perlogR=perlogR,
            ylog=ylog_num,
            fig_ax=[fig, axs[r, 1]],
        )
        pltdist.plot_distribs_in_select_range(
            time,
            pltdist.plot_domainmass_distribs,
            var_for_range,
            coordlims,
            superdrops,
            t2plts,
            volume,
            rspan,
            nbins,
            smoothsig=smoothsig,
            perlogR=perlogR,
            ylog=ylog_mass,
            fig_ax=[fig, axs[r, 2]],
        )
        for ax in axs[r, :]:
            ax.set_title(f"distribution between {coordlims[0]} < z < {coordlims[1]}")

    fig.tight_layout()
    if savename != "":
        fig.savefig(savename, dpi=400, bbox_inches="tight", facecolor="w", format="png")
        print("Figure .png saved as: " + str(savename))


# %%
### -------------------------------- MAIN ---------------------------------- ###
def main(
    path2CLEO,
    savefigpath,
    grid_filename,
    setupfile,
    dataset,
):
    import matplotlib.pyplot as plt

    from plotcleo import pltsds, pltmoms
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
    savename = savefigpath / "eurec4a1d_totnsupers.png"
    pltmoms.plot_totnsupers(time, totnsupers, savename=savename)
    plt.show()

    savename = savefigpath / "eurec4a1d_domainmassmoms.png"
    pltmoms.plot_domainmassmoments(time, massmoms, savename=savename)
    plt.show()

    nsample = 500
    savename = savefigpath / "eurec4a1d_randomsample.png"
    pltsds.plot_randomsample_superdrops(time, superdrops, nsample, savename=savename)
    plt.show()

    savename = savefigpath / "eurec4a1d_droplet_distributions.png"
    plot_droplet_distributions(config, gbxs, time, superdrops, savename=savename)
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
