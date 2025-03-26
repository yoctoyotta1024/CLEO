"""
----- CLEO -----
File: divfree2d_inputfiles.py
Project: divfreemotion
Created Date: Friday 17th November 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
Last Modified: Friday 3rd May 2024
Modified By: CB
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
Copyright (c) 2023 MPI-M, Clara Bayley
-----
File Description:
Script generates input files for divergence free motion of superdroplets in
a 2-D divergence free wind field.
"""

import sys


def main(
    path2CLEO,
    path2build,
    config_filename,
    grid_filename,
    initsupers_filename,
    thermofiles,
):
    import numpy as np

    sys.path.append(str(path2CLEO))  # for imports from pySD package

    from pySD import geninitconds
    from pySD.initsuperdropsbinary_src import (
        crdgens,
        rgens,
        dryrgens,
        probdists,
        attrsgen,
    )
    from pySD.thermobinary_src import thermogen, windsgen, thermodyngen

    ### ---------------------------------------------------------------- ###
    ### ----------------------- INPUT PARAMETERS ----------------------- ###
    ### ---------------------------------------------------------------- ###
    ### --- essential paths and filenames --- ###
    # path and filenames for creating initial SD conditions
    constants_filename = path2CLEO / "libs" / "cleoconstants.hpp"

    ### --- plotting initialisation figures --- ###
    isfigures = [
        True,
        True,
    ]  # booleans for [making, saving] initialisation figures # TODO(CB): move into args
    savefigpath = (
        path2build / "bin"
    )  # directory for saving figures # TODO(CB): move into args
    SDgbxs2plt = [
        0
    ]  # gbxindex of SDs to plot (nb. "all" can be very slow) # TODO(CB): move into args

    ### --- settings for 2-D gridbox boundaries --- ###
    zgrid = [0, 1500, 50]  # evenly spaced zhalf coords [zmin, zmax, zdelta] [m]
    xgrid = [0, 1500, 50]  # evenly spaced xhalf coords [m]
    ygrid = np.array([0, 20])  # array of yhalf coords [m]

    ### --- settings for initial superdroplets --- ###
    # settings for initial superdroplet coordinates
    zlim = 750  # max z coord of superdroplets
    npergbx = 8  # number of superdroplets per gridbox

    # [min, max] range of initial superdroplet radii (and implicitly solute masses)
    rspan = [3e-9, 3e-6]  # [m]

    # settings for initial superdroplet multiplicies
    # (from bimodal Lognormal distribution)
    geomeans = [0.02e-6, 0.15e-6]
    geosigs = [1.4, 1.6]
    scalefacs = [6e6, 4e6]
    numconc = np.sum(scalefacs)

    ### --- settings for 2D Thermodynamics --- ###
    PRESS = 100000  # [Pa]
    THETA = 298.15  # [K]
    qcond = 0.0  # [Kg/Kg]
    WMAX = 0.6  # [m/s]
    VVEL = None  # [m/s]
    Zlength = 1500  # [m]
    Xlength = 1500  # [m]
    qvapmethod = "sratio"
    Zbase = 750  # [m]
    sratios = [1.0, 1.0]  # s_ratio [below, above] Zbase
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

    ### ----- write thermodynamics binaries ----- ###
    thermog = thermogen.Simple2TierRelativeHumidity(
        config_filename,
        constants_filename,
        PRESS,
        THETA,
        qvapmethod,
        sratios,
        Zbase,
        qcond,
    )
    windsg = windsgen.Simple2DFlowField(
        config_filename, constants_filename, WMAX, Zlength, Xlength, VVEL
    )
    thermodyngen = thermodyngen.ThermodynamicsGenerator(thermog, windsg)
    geninitconds.generate_thermodynamics_conditions_fromfile(
        thermofiles,
        thermodyngen,
        config_filename,
        constants_filename,
        grid_filename,
        isfigures=isfigures,
        savefigpath=savefigpath,
    )

    ### ----- write initial superdroplets binary ----- ###
    nsupers = crdgens.nsupers_at_domain_base(
        grid_filename, constants_filename, npergbx, zlim
    )
    coord3gen = crdgens.SampleCoordGen(True)  # sample coord3 randomly
    coord1gen = crdgens.SampleCoordGen(True)  # sample coord1 randomly
    coord2gen = None  # do not generate superdroplet coord2s
    xiprobdist = probdists.LnNormal(geomeans, geosigs, scalefacs)
    radiigen = rgens.SampleLog10RadiiGen(rspan)  # randomly sample radii from rspan [m]
    dryradiigen = dryrgens.ScaledRadiiGen(1.0)

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
    ### args = path2CLEO, path2build, config_filename, grid_filename, initsupers_filename, thermofiles
    main(*sys.argv[1:])
