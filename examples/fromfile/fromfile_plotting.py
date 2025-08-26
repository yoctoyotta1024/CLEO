"""
Copyright (c) 2024 MPI-M, Clara Bayley


----- CLEO -----
File: fromfile_plotting.py
Project: fromfile
Created Date: Thursday 30th January 2025
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
Script plots results of 3D example with time varying thermodynamics
read from binary files.
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
        default="/home/m/m300950/CLEO/build_fromfile/bin/ntasks4",
    )
    parser.add_argument(
        "--grid_filename",
        type=Path,
        help="Absolute path to gridbox boundaries file",
        default="/home/m/m300950/CLEO/build_fromfile/share/fromfile_dimlessGBxboundaries.dat",
    )
    parser.add_argument(
        "--setupfile",
        type=Path,
        help="Absolute path to setup file",
        default="/home/m/m300950/CLEO/build_fromfile/bin/ntasks4/fromfile_setup.txt",
    )
    parser.add_argument(
        "--dataset",
        type=Path,
        help="Absolute path to dataset",
        default="/home/m/m300950/CLEO/build_fromfile/bin/ntasks4/fromfile_sol.zarr",
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

    from src import plot_output_thermo
    from plotssrc import pltsds, pltmoms
    from cleopy.sdmout_src import pyzarr, pysetuptxt, pygbxsdat

    # read in constants and intial setup from setup .txt file
    config = pysetuptxt.get_config(setupfile, nattrs=3, isprint=True)
    consts = pysetuptxt.get_consts(setupfile, isprint=True)
    gbxs = pygbxsdat.get_gridboxes(grid_filename, consts["COORD0"], isprint=True)

    time = pyzarr.get_time(dataset)
    superdrops = pyzarr.get_supers(dataset, consts)
    maxnsupers = pyzarr.get_totnsupers(dataset)
    thermo, winds = pyzarr.get_thermodata(
        dataset, config["ntime"], gbxs["ndims"], consts, getwinds=True
    )

    # plot super-droplet results
    savename = savefigpath / "fromfile_maxnsupers_validation.png"
    pltmoms.plot_totnsupers(time, maxnsupers, savename=savename)
    plt.show()

    nsample = 1000
    savename = savefigpath / "fromfile_motion2d_validation.png"
    pltsds.plot_randomsample_superdrops_2dmotion(
        superdrops,
        nsample,
        savename=savename,
        arrows=False,
        israndom=False,
    )
    plt.show()

    # plot thermodynamics results
    plot_output_thermo.plot_domain_thermodynamics_timeseries(
        time,
        gbxs,
        thermo,
        winds,
        savedir=savefigpath,
        do_show=True,
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
