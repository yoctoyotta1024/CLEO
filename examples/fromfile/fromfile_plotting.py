"""
Copyright (c) 2024 MPI-M, Clara Bayley


----- CLEO -----
File: fromfile_plotting.py
Project: fromfile
Created Date: Thursday 30th January 2025
Author: Clara Bayley (CB)
Additional Contributors:
-----
Last Modified: Thursday 30th January 2025
Modified By: CB
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
Script plots results of 3D example with time varying thermodynamics
read from binary files.
"""

import sys


def main(
    path2CLEO,
    grid_filename,
    setupfile,
    dataset,
    savefigpath,
):
    sys.path.append(str(path2CLEO))  # imports from pySD
    sys.path.append(
        str(path2CLEO / "examples" / "exampleplotting")
    )  # imports from example plots package

    from src import plot_output_thermo
    from plotssrc import pltsds, pltmoms
    from pySD.sdmout_src import pyzarr, pysetuptxt, pygbxsdat

    # read in constants and intial setup from setup .txt file
    config = pysetuptxt.get_config(setupfile, nattrs=3, isprint=True)
    consts = pysetuptxt.get_consts(setupfile, isprint=True)
    gbxs = pygbxsdat.get_gridboxes(grid_filename, consts["COORD0"], isprint=True)

    time = pyzarr.get_time(dataset)
    sddata = pyzarr.get_supers(dataset, consts)
    maxnsupers = pyzarr.get_totnsupers(dataset)
    thermo, winds = pyzarr.get_thermodata(
        dataset, config["ntime"], gbxs["ndims"], consts, getwinds=True
    )

    # plot super-droplet results
    savename = savefigpath / "fromfile_maxnsupers_validation.png"
    pltmoms.plot_totnsupers(time, maxnsupers, savename=savename)

    nsample = 1000
    savename = savefigpath / "fromfile_motion2d_validation.png"
    pltsds.plot_randomsample_superdrops_2dmotion(
        sddata,
        config["maxnsupers"],
        nsample,
        savename=savename,
        arrows=False,
        israndom=False,
    )

    # plot thermodynamics results
    plot_output_thermo.plot_domain_thermodynamics_timeseries(
        time, gbxs, thermo, winds, savedir=savefigpath
    )


if __name__ == "__main__":
    ### args = path2CLEO, grid_filename, setupfile, dataset, savefigpath,
    main(*sys.argv[1:])
