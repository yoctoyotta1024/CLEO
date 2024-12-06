"""
Copyright (c) 2024 MPI-M, Clara Bayley


----- CLEO -----
File: bubble3d_inputfiles.py
Project: bubble3d
Created Date: Friday 17th November 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
Last Modified: Tuesday 9th July 2024
Modified By: CB
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
Script generates input files for 3D example with time varying
thermodynamics read from ICON output of bubble test case by YAC
"""

import sys
from pathlib import Path


def get_zgrid(icon_grid_file, num_vertical_levels):
    """returns zgrid for CLEO gridfile with same vertical levels as ICON grid file"""
    import numpy as np
    import xarray as xr

    grid = xr.open_dataset(icon_grid_file)
    idx2 = int(grid.height.values[-1])
    idx1 = int(idx2 - num_vertical_levels - 1)
    zhalf = grid.zghalf.values[idx1:idx2, 0]  # [m]
    zgrid = np.flip(zhalf)

    return zgrid  # [m]


def main(
    path2CLEO,
    path2build,
    config_filename,
    grid_filename,
    initsupers_filename,
    icon_grid_file,
    SDgbxs2plt,
):
    sys.path.append(str(path2CLEO))  # for imports from pySD package

    from pySD import geninitconds
    from pySD.initsuperdropsbinary_src import (
        crdgens,
        rgens,
        dryrgens,
        probdists,
        attrsgen,
    )

    ### ---------------------------------------------------------------- ###
    ### ----------------------- INPUT PARAMETERS ----------------------- ###
    ### ---------------------------------------------------------------- ###
    ### --- essential paths and filenames --- ###
    # path and filenames for creating initial SD conditions
    constants_filename = path2CLEO / Path("libs/cleoconstants.hpp")

    ### --- plotting initialisation figures --- ###
    # booleans for [making, saving] initialisation figures
    isfigures = [True, True]  # TODO(CB): move into args
    savefigpath = path2build / Path("bin")  # binpath # TODO(CB): move into args

    ### --- settings for 3-D gridbox boundaries --- ###
    num_vertical_levels = 24  # TODO(CB): move to config file (?)
    zgrid = get_zgrid(icon_grid_file, num_vertical_levels)  # [m]
    xgrid = [
        0,
        30000,
        2500,
    ]  # evenly spaced xhalf coords [m] # distance must match longitude in config file
    ygrid = [
        0,
        6250,
        1250,
    ]  # evenly spaced xhalf coords [m] # distance must match latitudes in config file

    ### --- settings for initial superdroplets --- ###
    # settings for initial superdroplet coordinates
    zlim = 1000  # max z coord of superdroplets
    npergbx = 2  # number of superdroplets per gridbox

    monor = 1e-6  # all SDs have this same radius [m]
    dryr_sf = 1.0  # scale factor for dry radii [m]
    numconc = 5e8  # total no. conc of real droplets [m^-3]
    randcoord = False  # sample SD spatial coordinates randomly or not
    ### ---------------------------------------------------------------- ###
    ### ---------------------------------------------------------------- ###

    if path2CLEO == path2build:
        raise ValueError("build directory cannot be CLEO")

    ### ---------------------------------------------------------------- ###
    ### ------------------- BINARY FILES GENERATION--------------------- ###
    ### ---------------------------------------------------------------- ###
    ### ----- write gridbox boundaries binary ----- ###
    geninitconds.generate_gridbox_boundaries(
        grid_filename,
        zgrid,
        xgrid,
        ygrid,
        constants_filename,
        isfigures=isfigures,
        savefigpath=savefigpath,
    )

    ### ----- write initial superdroplets binary ----- ###
    nsupers = crdgens.nsupers_at_domain_base(
        grid_filename, constants_filename, npergbx, zlim
    )
    radiigen = rgens.MonoAttrGen(monor)  # all SDs have the same radius [m]
    dryradiigen = dryrgens.ScaledRadiiGen(dryr_sf)  # dryradii are 1/sf of radii [m]
    coord3gen = crdgens.SampleCoordGen(randcoord)  # (not) random coord3 of SDs
    coord1gen = crdgens.SampleCoordGen(randcoord)  # (not) random coord1 of SDs
    coord2gen = crdgens.SampleCoordGen(randcoord)  # (not) random coord2 of SDs
    xiprobdist = probdists.DiracDelta(monor)  # monodisperse droplet probability distrib

    initattrsgen = attrsgen.AttrsGenerator(
        radiigen, dryradiigen, xiprobdist, coord3gen, coord1gen, coord2gen
    )
    geninitconds.generate_initial_superdroplet_conditions(
        initattrsgen,
        initsupers_filename,
        config_filename,
        constants_filename,
        grid_filename,
        nsupers,
        numconc,
        isfigures=isfigures,
        savefigpath=savefigpath,
        gbxs2plt=SDgbxs2plt,
    )
    ### ---------------------------------------------------------------- ###
    ### ---------------------------------------------------------------- ###


if __name__ == "__main__":
    ### args = path2CLEO, path2build, config_filename, grid_filename, initsupers_file, icon_grid_file, SDgbxs2plt
    main(*sys.argv[1:])
