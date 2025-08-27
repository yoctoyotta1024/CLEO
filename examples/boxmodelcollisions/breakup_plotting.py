"""
Copyright (c) 2025 MPI-M, Clara Bayley


----- CLEO -----
File: breakup_plotting.py
Project: boxmodelcollisions
Created Date: Friday 22nd August 2025
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
Script for plotting results of CLEO 0-D box model for collisions using
selected collision kernels with breakup (e.g. Low and Lists's) in a way
comparable to Shima et al. 2009 Fig. 2(b)
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
        default="/home/m/m300950/CLEO/build_colls0d/breakup/bin",
    )
    parser.add_argument(
        "--grid_filename",
        type=Path,
        help="Absolute path to gridbox boundaries file",
        default="/home/m/m300950/CLEO/build_colls0d/breakup/share/breakup_dimlessGBxboundaries.dat",
    )
    parser.add_argument(
        "--setupfiles",
        nargs="*",
        type=Path,
        help="Absolute path to setup files",
        default=[
            "/home/m/m300950/CLEO/build_colls0d/breakup/bin/breakup_long_setup.txt",
        ],
    )
    parser.add_argument(
        "--datasets",
        nargs="*",
        type=Path,
        help="Absolute path to dataset",
        default=[
            "/home/m/m300950/CLEO/build_colls0d/breakup/bin/breakup_long_sol.zarr"
        ],
    )
    parser.add_argument(
        "--kernels",
        nargs="*",
        choices=["long", "lowlist", "szakallurbich", "testikstraub"],
        type=str,
        help="kernel examples to plot",
        default=["long"],
    )
    return parser.parse_args()


def get_results(path2CLEO, setupfile, dataset):
    """read in time, superdrops and setup for given dataset and setupfile"""
    from cleopy.sdmout_src import pyzarr, pysetuptxt

    # read in constants and intial setup from setup .txt file
    config = pysetuptxt.get_config(setupfile, nattrs=3, isprint=True)
    consts = pysetuptxt.get_consts(setupfile, isprint=True)

    # read in seom data from dataset
    time = pyzarr.get_time(dataset)
    superdrops = pyzarr.get_supers(dataset, consts)

    return config, consts, time, superdrops


def plot_onekernel_results(
    path2CLEO,
    grid_filename,
    setupfile,
    dataset,
    xlims,
    t2plts,
    savename,
):
    import sys
    import matplotlib.pyplot as plt

    sys.path.append(
        str(path2CLEO / "examples" / "exampleplotting")
    )  # imports from example plots package
    from plotssrc import shima2009fig
    from cleopy.sdmout_src import pygbxsdat

    # read in data
    config, consts, time, superdrops = get_results(path2CLEO, setupfile, dataset)
    gbxs = pygbxsdat.get_gridboxes(grid_filename, consts["COORD0"], isprint=True)

    # make and save plot
    numconc = 2**23  # total no. conc of real droplets [m^-3]
    volexpr0 = 30.531e-6  # peak of volume exponential distribution [m]
    smoothsigconst = 0.62
    smoothsig = smoothsigconst * (
        config["maxnsupers"] ** (-1 / 5)
    )  # = ~0.2 for guassian smoothing
    withgol, plotwitherr = False, False

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
    path2CLEO, grid_filename, kernels, setupfiles, datasets, xlims, t2plts, savename
):
    import sys
    import awkward as ak
    import matplotlib.pyplot as plt

    sys.path.append(
        str(path2CLEO / "examples" / "exampleplotting")
    )  # imports from example plots package
    from plotssrc import shima2009fig
    from cleopy.sdmout_src import pygbxsdat

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
    for ker, set, dat in zip(kernels, setupfiles, datasets):
        # read in data
        config, consts, time, superdrops = get_results(path2CLEO, set, dat)
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
            if ker == "long":
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
                    label=styles[ker][0],
                    color=styles[ker][1],
                    linestyle=linestyles[n],
                    linewidth=styles[ker][2],
                )
            else:
                ax.plot(
                    hcens,
                    hist,
                    color=styles[ker][1],
                    linestyle=linestyles[n],
                    linewidth=styles[ker][2],
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


def main(path2CLEO, savefigpath, grid_filename, setupfiles, datasets, kernels):
    assert len(setupfiles) == len(datasets) and len(kernels) == len(
        datasets
    ), "number of items in lists of setupfiles, datasets and kernels should be equal"

    for ker, set, dat in zip(kernels, setupfiles, datasets):
        t2plts = [0, 600, 1200, 1800, 2400]
        xlims = [10, 5000]
        savename = savefigpath / f"{ker}_validation.png"
        plot_onekernel_results(
            path2CLEO,
            grid_filename,
            set,
            dat,
            xlims,
            t2plts,
            savename,
        )

    t2plts = [0, 600, 1200, 1800, 2400]
    xlims = [10, 5000]
    savename = savefigpath / "breakup_validation.png"
    plot_allkernels_results(
        path2CLEO, grid_filename, kernels, setupfiles, datasets, xlims, t2plts, savename
    )


# %%
### --------------------------- RUN PROGRAM -------------------------------- ###
if __name__ == "__main__":
    args = parse_arguments()
    main(
        args.path2CLEO,
        args.savefigpath,
        args.grid_filename,
        args.setupfiles,
        args.datasets,
        args.kernels,
    )
