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


def main(path2CLEO, path2build, configfile, gridfile, initSDsfile, thermofiles):
    import numpy as np
    import matplotlib.pyplot as plt

    sys.path.append(str(path2CLEO))  # for imports from pySD package

    from pySD.gbxboundariesbinary_src import read_gbxboundaries as rgrid
    from pySD.gbxboundariesbinary_src import create_gbxboundaries as cgrid
    from pySD.initsuperdropsbinary_src import (
        crdgens,
        rgens,
        dryrgens,
        probdists,
        attrsgen,
    )
    from pySD.initsuperdropsbinary_src import create_initsuperdrops as csupers
    from pySD.initsuperdropsbinary_src import read_initsuperdrops as rsupers
    from pySD.thermobinary_src import thermogen
    from pySD.thermobinary_src import create_thermodynamics as cthermo
    from pySD.thermobinary_src import read_thermodynamics as rthermo

    ### ---------------------------------------------------------------- ###
    ### ----------------------- INPUT PARAMETERS ----------------------- ###
    ### ---------------------------------------------------------------- ###
    ### --- essential paths and filenames --- ###
    # path and filenames for creating initial SD conditions
    constsfile = path2CLEO / "libs" / "cleoconstants.hpp"

    ### --- plotting initialisation figures --- ###
    isfigures = [True, True]  # booleans for [making, saving] initialisation figures
    savefigpath = path2build / "bin"  # directory for saving figures
    SDgbxs2plt = [0]  # gbxindex of SDs to plot (nb. "all" can be very slow)

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
    PRESS0 = 100000  # [Pa]
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
    cgrid.write_gridboxboundaries_binary(gridfile, zgrid, xgrid, ygrid, constsfile)
    rgrid.print_domain_info(constsfile, gridfile)

    ### ----- write thermodynamics binaries ----- ###
    thermodyngen = thermogen.SimpleThermo2DFlowField(
        configfile,
        constsfile,
        PRESS0,
        THETA,
        qvapmethod,
        sratios,
        Zbase,
        qcond,
        WMAX,
        Zlength,
        Xlength,
        VVEL,
    )
    cthermo.write_thermodynamics_binary(
        thermofiles, thermodyngen, configfile, constsfile, gridfile
    )

    ### ----- write initial superdroplets binary ----- ###
    nsupers = crdgens.nsupers_at_domain_base(gridfile, constsfile, npergbx, zlim)
    coord3gen = crdgens.SampleCoordGen(True)  # sample coord3 randomly
    coord1gen = crdgens.SampleCoordGen(True)  # sample coord1 randomly
    coord2gen = None  # do not generate superdroplet coord2s
    xiprobdist = probdists.LnNormal(geomeans, geosigs, scalefacs)
    radiigen = rgens.SampleLog10RadiiGen(rspan)  # randomly sample radii from rspan [m]
    dryradiigen = dryrgens.ScaledRadiiGen(1.0)

    initattrsgen = attrsgen.AttrsGenerator(
        radiigen, dryradiigen, xiprobdist, coord3gen, coord1gen, coord2gen
    )
    csupers.write_initsuperdrops_binary(
        initSDsfile, initattrsgen, configfile, constsfile, gridfile, nsupers, numconc
    )

    ### ----- show (and save) plots of binary file data ----- ###
    if isfigures[0]:
        rgrid.plot_gridboxboundaries(constsfile, gridfile, savefigpath, isfigures[1])
        rthermo.plot_thermodynamics(
            constsfile, configfile, gridfile, thermofiles, savefigpath, isfigures[1]
        )
        rsupers.plot_initGBxs_distribs(
            configfile,
            constsfile,
            initSDsfile,
            gridfile,
            savefigpath,
            isfigures[1],
            SDgbxs2plt,
        )
        plt.close()
    ### ---------------------------------------------------------------- ###
    ### ---------------------------------------------------------------- ###


if __name__ == "__main__":
    ### args = path2CLEO, path2build, configfile, binpath, gridfile, initSDsfile, thermofiles
    main(*sys.argv[1:])
