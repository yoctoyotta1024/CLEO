"""
Copyright (c) 2025 MPI-M, Clara Bayley


----- CLEO -----
File: shima2009_plotting.py
Project: boxmodelcollisions
Created Date: Friday 22nd August 2025
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
Script for plotting results of CLEO 0-D box model for condensation and collisions
using the Long kernel in a somewhat comparable way to Shima et al. 2009 Fig. 2
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
        default="/home/m/m300950/CLEO/build_colls0d/shima2009/bin",
    )
    parser.add_argument(
        "--grid_filename",
        type=Path,
        help="Absolute path to gridbox boundaries file",
        default="/home/m/m300950/CLEO/build_colls0d/shima2009/share/shima2009_dimlessGBxboundaries.dat",
    )
    parser.add_argument(
        "--setupfile",
        type=Path,
        help="Absolute path to setup file",
        default="/home/m/m300950/CLEO/build_colls0d/shima2009/bin/shima2009_golovin_setup.txt",
    )
    parser.add_argument(
        "--dataset",
        type=Path,
        help="Absolute path to dataset",
        default="/home/m/m300950/CLEO/build_colls0d/shima2009/bin/shima2009_golovin_sol.zarr/",
    )
    return parser.parse_args()


# %%
def plot_droplet_distributions(
    config, gbxs, time, superdrops, t2plts=None, savename=""
):
    import awkward as ak
    import matplotlib.pyplot as plt
    from plotcleo import pltdist

    fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(16, 8), sharex=True)

    if t2plts is None:
        t2plts = [time.secs[0], time.secs[10], time.secs[-1]]

    volume = gbxs["domainvol"]
    rspan = [ak.min(superdrops["radius"]) * 0.9, ak.max(superdrops["radius"]) * 1.1]
    nbins = 100
    smoothsig = False
    perlogR = False
    ylog_nsupers = False
    ylog_num = False
    ylog_mass = True

    fig, ax = pltdist.plot_domainnsupers_distribs(
        time,
        superdrops,
        t2plts,
        volume,
        rspan,
        nbins,
        smoothsig=smoothsig,
        perlogR=perlogR,
        ylog=ylog_nsupers,
        fig_ax=[fig, axs[0]],
        savename=savename,
    )

    fig, ax = pltdist.plot_domainnumconc_distribs(
        time,
        superdrops,
        t2plts,
        volume,
        rspan,
        nbins,
        smoothsig=smoothsig,
        perlogR=perlogR,
        ylog=ylog_num,
        fig_ax=[fig, axs[1]],
        savename=savename,
    )

    fig, ax = pltdist.plot_domainmass_distribs(
        time,
        superdrops,
        t2plts,
        volume,
        rspan,
        nbins,
        smoothsig=smoothsig,
        perlogR=perlogR,
        ylog=ylog_mass,
        fig_ax=[fig, axs[2]],
        savename=savename,
    )

    fig.tight_layout()
    if savename != "":
        fig.savefig(savename, dpi=400, bbox_inches="tight", facecolor="w", format="png")
        print("Figure .png saved as: " + str(savename))


def main(path2CLEO, savefigpath, grid_filename, setupfile, dataset):
    import matplotlib.pyplot as plt

    from plotcleo import pltmoms
    from cleopy.sdmout_src import pyzarr, pysetuptxt, pygbxsdat

    # plot settings
    t2plts = [0, 200, 600, 1200, 1800]

    # read in constants and intial setup from setup .txt file
    config = pysetuptxt.get_config(setupfile, nattrs=3, isprint=True)
    consts = pysetuptxt.get_consts(setupfile, isprint=True)
    gbxs = pygbxsdat.get_gridboxes(grid_filename, consts["COORD0"], isprint=True)

    time = pyzarr.get_time(dataset)
    superdrops = pyzarr.get_supers(dataset, consts)
    massmoms = pyzarr.get_massmoms(dataset, config["ntime"], gbxs["ndims"])

    savename = savefigpath / "boxmodel_condcoll_domainmassmoms.png"
    pltmoms.plot_domainmassmoments(time, massmoms, savename=savename)
    plt.show()

    savename = savefigpath / "boxmodel_condcoll_dsd_evolution.png"
    plot_droplet_distributions(
        config, gbxs, time, superdrops, t2plts=t2plts, savename=savename
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
