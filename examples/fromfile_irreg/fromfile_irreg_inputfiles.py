"""
Copyright (c) 2024 MPI-M, Clara Bayley


----- CLEO -----
File: fromfile_irreg_inputfiles.py
Project: fromfile_irreg
Created Date: Wednesday 11th September 2024
Author: Clara Bayley (CB)
Additional Contributors:
-----
Last Modified: Wednesday 11th September 2024
Modified By: CB
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
Script generates input files for 3D example with irregular grid and
time varying thermodynamics read from binary files.
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

    from src import gen_input_thermo
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
    constants_filename = path2CLEO / "libs" / "cleoconstants.hpp"

    ### --- plotting initialisation figures --- ###
    # booleans for [making, saving] initialisation figures
    isfigures = [True, True]  # TODO(CB): move into args
    savefigpath = path2build / "bin"  # binpath # TODO(CB): move into args

    ### --- settings for 3-D irregular gridbox boundaries --- ###
    zgrid = np.array([0, 20, 30, 45, 60, 80, 90, 120, 140, 180, 360, 500, 1000, 1500])
    xgrid = np.array(
        [0, 33, 205, 440, 650, 915, 1033, 1100, 1300, 1450, 1500]
    )  # evenly spaced xhalf coords [m]
    ygrid = np.array([0, 10, 75, 100, 150, 200, 300])  # array of yhalf coords [m]

    ### --- settings for initial superdroplets --- ###
    # settings for initial superdroplet coordinates
    zlim = 1000  # max z coord of superdroplets
    npergbx = 2  # number of superdroplets per gridbox

    monor = 1e-6  # all SDs have this same radius [m]
    dryr_sf = 1.0  # scale factor for dry radii [m]
    numconc = 5e8  # total no. conc of real droplets [m^-3]
    randcoord = False  # sample SD spatial coordinates randomly or not

    ### --- settings for 2D Thermodynamics --- ###
    PRESSz0 = 101500  # [Pa]
    TEMPz0 = 300  # [K]
    qvapz0 = 0.05  # [Kg/Kg]
    qcondz0 = 0.001  # [Kg/Kg]
    WMAX = 1.5  # [m/s]
    Zlength = 1500  # [m]
    Xlength = 1500  # [m]
    VMAX = 1.0  # [m/s]
    Ylength = 300  # [m]
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
    thermodyngen = gen_input_thermo.TimeVarying3DThermo(
        PRESSz0, TEMPz0, qvapz0, qcondz0, WMAX, Zlength, Xlength, VMAX, Ylength
    )
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
    )
    ### ---------------------------------------------------------------- ###
    ### ---------------------------------------------------------------- ###


if __name__ == "__main__":
    ### args = path2CLEO, path2build, config_filename, grid_filename, initsupers_filename, thermofiles
    main(*sys.argv[1:])
