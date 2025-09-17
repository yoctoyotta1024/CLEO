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
