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
Script for plotting results of CLEO 0-D box model for collisions using the
Golovin or Long (hydrodynamic) kernel in a way comparable to Shima et al. 2009 Fig. 2
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
    parser.add_argument(
        "--kernel",
        type=str,
        choices=["golovin", "long1", "long2"],
        help="name of kernel example to plot",
        default="golovin",
    )
    return parser.parse_args()


def plot_shima2009_plot_validation_figure(
    path2CLEO,
    grid_filename,
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
    import sys
    import matplotlib.pyplot as plt

    sys.path.append(
        str(path2CLEO / "examples" / "exampleplotting")
    )  # imports from example plots package
    from plotssrc import shima2009fig
    from cleopy.sdmout_src import pyzarr, pysetuptxt, pygbxsdat

    # read in constants and intial setup from setup .txt file
    config = pysetuptxt.get_config(setupfile, nattrs=3, isprint=True)
    consts = pysetuptxt.get_consts(setupfile, isprint=True)
    gbxs = pygbxsdat.get_gridboxes(grid_filename, consts["COORD0"], isprint=True)

    time = pyzarr.get_time(dataset)
    superdrops = pyzarr.get_supers(dataset, consts)

    # make and save plot
    smoothsig = smoothsigconst * (
        config["maxnsupers"] ** (-1 / 5)
    )  # = ~0.2 for guassian smoothing
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


def main(path2CLEO, savefigpath, grid_filename, setupfile, dataset, kernel):
    if "golovin" == kernel or "long1" == kernel:
        numconc = 2**23  # total no. conc of real droplets [m^-3]
        volexpr0 = 30.531e-6  # peak of volume exponential distribution [m]
        smoothsigconst = 0.62
        xlims = [10, 5000]
        if "golovin" == kernel:
            t2plts = [0, 1200, 2400, 3600]
            withgol, plotwitherr = True, True
            savename = savefigpath / "golovin_validation.png"
        elif "long1" == kernel:
            t2plts = [0, 600, 1200, 1800]
            withgol, plotwitherr = False, False
            savename = savefigpath / "long_validation_1.png"
    elif "long2" == kernel:
        numconc = 27 * 2**23  # total no. conc of real droplets [m^-3]
        volexpr0 = 10.117e-6  # peak of volume exponential distribution [m]
        smoothsigconst = 1.0
        xlims = [1, 5000]
        t2plts = [0, 1200, 1800, 2400, 3600]
        withgol, plotwitherr = False, False
        savename = savefigpath / "long_validation_2.png"
    else:
        raise ValueError(
            "kernel for examples not recognised, please choose from: golovin, long1 or long2"
        )

    plot_shima2009_plot_validation_figure(
        path2CLEO,
        grid_filename,
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
        args.kernel,
    )
