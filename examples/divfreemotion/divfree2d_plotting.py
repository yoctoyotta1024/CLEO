"""
Copyright (c) 2025 MPI-M, Clara Bayley


----- CLEO -----
File: divfree2d_plotting.py
Project: divfreemotion
Created Date: Friday 22nd August 2025
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
Script plots motion of superdroplets in 2-D divergence free flow example
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
        default="/home/m/m300950/CLEO/build_divfree2d/bin",
    )
    parser.add_argument(
        "--grid_filename",
        type=Path,
        help="Absolute path to gridbox boundaries file",
        default="/home/m/m300950/CLEO/build_divfree2d/share/divfree2d_dimlessGBxboundaries.dat",
    )
    parser.add_argument(
        "--setupfile",
        type=Path,
        help="Absolute path to setup file",
        default="/home/m/m300950/CLEO/build_divfree2d/bin/divfree2d_setup.txt",
    )
    parser.add_argument(
        "--dataset",
        type=Path,
        help="Absolute path to dataset",
        default="/home/m/m300950/CLEO/build_divfree2d/bin/divfree2d_sol.zarr",
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
    import matplotlib.pyplot as plt

    sys.path.append(
        str(path2CLEO / "examples" / "exampleplotting")
    )  # imports from example plots package

    from plotssrc import pltsds, pltmoms
    from cleopy.sdmout_src import pyzarr, pysetuptxt

    # read in constants and dataset
    consts = pysetuptxt.get_consts(setupfile, isprint=True)
    time = pyzarr.get_time(dataset)
    superdrops = pyzarr.get_supers(dataset, consts)
    maxnsupers = pyzarr.get_totnsupers(dataset)

    # 4. plot results
    savename = savefigpath / "divfree2d_maxnsupers_validation.png"
    pltmoms.plot_totnsupers(time, maxnsupers, savename=savename)
    plt.show()

    nsample = 500
    savename = savefigpath / "divfree2d_motion2d_validation.png"
    pltsds.plot_randomsample_superdrops_2dmotion(
        superdrops, nsample, savename=savename, arrows=False
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
