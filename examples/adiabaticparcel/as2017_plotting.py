"""
Copyright (c) 2025 MPI-M, Clara Bayley


----- CLEO -----
File: as2017_plotting.py
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
Script plots results for adiabatic parcel example similar to Figure 5 of "On
the CCN (de)activation nonlinearities" S. Arabas and S. Shima 2017 to show
example of adaibatic parcel expansion and contraction.
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
        default="/home/m/m300950/CLEO/build_adia0d/as2017/bin",
    )
    parser.add_argument(
        "--grid_filename",
        type=Path,
        help="Absolute path to gridbox boundaries file",
        default="/home/m/m300950/CLEO/build_adia0d/as2017/share/as2017_dimlessGBxboundaries.dat",
    )
    parser.add_argument(
        "--setupfiles",
        nargs="*",
        type=Path,
        help="Absolute path to setup files",
        default=[
            "/home/m/m300950/CLEO/build_adia0d/as2017/bin/as2017_setup_run0.txt",
        ],
    )
    parser.add_argument(
        "--datasets",
        nargs="*",
        type=Path,
        help="Absolute path to dataset",
        default=["/home/m/m300950/CLEO/build_adia0d/as2017/bin/as2017_sol_run0.zarr"],
    )
    parser.add_argument(
        "--runnums",
        nargs="*",
        type=int,
        choices=[0, 1, 2, 3, 4, 5, 6, 7, 8],
        help="Run number of each setupfile, dataset pair",
        default=[0],
    )
    return parser.parse_args()


# %%
### ------------------------- FUNCTION DEFINITIONS ------------------------- ###
def displacement(time, w_avg, thalf):
    """displacement z given velocity, w, is sinusoidal
    profile: w = w_avg * pi/2 * np.sin(np.pi * t/thalf)
    where wmax = pi/2*w_avg and tauhalf = thalf/pi."""
    import numpy as np

    zmax = w_avg / 2 * thalf
    z = zmax * (1 - np.cos(np.pi * time / thalf))
    return z


def get_results(path2CLEO, grid_filename, setupfile, dataset):
    """read in required data from given dataset and setupfile"""
    from cleopy.sdmout_src import pyzarr, pysetuptxt, pygbxsdat

    config = pysetuptxt.get_config(setupfile, nattrs=3, isprint=True)
    consts = pysetuptxt.get_consts(setupfile, isprint=True)
    gbxs = pygbxsdat.get_gridboxes(grid_filename, consts["COORD0"], isprint=True)

    # read in output Xarray data
    time = pyzarr.get_time(dataset).secs
    thermo = pyzarr.get_thermodata(dataset, len(time), gbxs["ndims"], consts)
    superdrops = pyzarr.get_supers(dataset, consts)
    zprof = displacement(time, config["W_avg"], config["TAU_half"])

    return config, gbxs, time, thermo, superdrops, zprof


# %%
### -------------------------------- MAIN ---------------------------------- ###
def main(path2CLEO, savefigpath, grid_filename, setupfiles, datasets, runnums):
    import sys
    import numpy as np
    import matplotlib.pyplot as plt

    sys.path.append(
        str(path2CLEO / "examples" / "exampleplotting")
    )  # imports from example plots package
    from plotssrc import as2017fig

    assert (
        len(setupfiles) == len(runnums) and len(datasets) == len(runnums)
    ), "number of items in lists of setupfiles and datasets should equal number of runnums"

    lwdths = {  # config_params number: linewidth
        0: 2,
        1: 1,
        2: 0.5,
    }

    fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(12, 16))
    r_initial_prev, numconc_prev = 0, 0
    for r, runnum in enumerate(runnums):
        icond = (
            runnum // 3
        )  # initial conditions number (3 = number of different params dicts)
        pnum = (
            runnum % 3
        )  # config_params number (3 = number of different initial conditions)
        axs = axes[:, icond]
        lwdth = lwdths[pnum]

        ### load results
        config, gbxs, time, thermo, superdrops, zprof = get_results(
            path2CLEO, grid_filename, setupfiles[r], datasets[r]
        )
        supersat = thermo.supersaturation()
        numconc = np.sum(superdrops["xi"][0]) / gbxs["domainvol"] / 1e6  # [/cm^3]

        attrs = ["radius", "xi", "msol"]
        sd0 = superdrops.sample("sdId", sample_values=0, variables2sample=attrs)
        radius = sd0["radius"][0]
        msol_initial = sd0["msol"][0][0]

        ### plot results
        wlab = "<w> = {:.1f}".format(config["W_avg"] * 100) + "cm s$^{-1}$"
        as2017fig.condensation_validation_subplots(
            axs, time, radius, supersat[:, 0, 0, 0], zprof, lwdth=lwdth, lab=wlab
        )

        ### save figure
        as2017fig.plot_kohlercurve_with_criticalpoints(
            axs[1],
            radius,
            msol_initial,
            thermo.temp[0, 0, 0, 0],
            superdrops.IONIC(),
            superdrops.MR_SOL(),
        )

        if pnum == 0:
            r_initial = radius[0]
            textlab = (
                "N = "
                + str(numconc)
                + "cm$^{-3}$\n"
                + "r$_{dry}$ = "
                + "{:.2g}\u03BCm\n".format(r_initial)
            )
            axs[0].text(0.03, 0.875, textlab, transform=axs[0].transAxes)
            axs[1].legend(loc="upper left")
            numconc_prev = numconc
            r_initial_prev = r_initial
        else:
            assert (
                numconc == numconc_prev and r_initial == r_initial_prev
            ), "all plots on same column of figure should have same numconc and initial radius"

    for ax in axes[0, :]:
        ax.legend(loc="lower right", fontsize=10)
        ax.set_xlim([-1, 1])
        ax.set_ylim([0, 150])
    for ax in axes[1:, :].flatten():
        ax.set_xlim([0.125, 10])
        ax.set_xscale("log")
    for ax in axes[1, :]:
        ax.set_ylim([-1, 1])
    for ax in axes[2, :]:
        ax.set_ylim([5, 75])

    fig.tight_layout()

    savename = savefigpath / "as2017fig.png"
    fig.savefig(savename, dpi=400, bbox_inches="tight", facecolor="w", format="png")
    print("Figure .png saved as: " + str(savename))
    plt.show()


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
        args.runnums,
    )
