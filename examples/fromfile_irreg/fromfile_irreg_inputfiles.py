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


def main(path2CLEO, path2build, configfile, gridfile, initSDsfile, thermofile):
    import numpy as np
    import matplotlib.pyplot as plt

    sys.path.append(path2CLEO)  # for imports from pySD package

    from src import gen_input_thermo
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
    from pySD.thermobinary_src import create_thermodynamics as cthermo
    from pySD.thermobinary_src import read_thermodynamics as rthermo

    ### ---------------------------------------------------------------- ###
    ### ----------------------- INPUT PARAMETERS ----------------------- ###
    ### ---------------------------------------------------------------- ###
    ### --- essential paths and filenames --- ###
    # path and filenames for creating initial SD conditions
    constsfile = path2CLEO + "/libs/cleoconstants.hpp"

    ### --- plotting initialisation figures --- ###
    # booleans for [making, saving] initialisation figures
    isfigures = [True, True]
    savefigpath = path2build + "/bin/"  # directory for saving figures

    ### --- settings for 2-D gridbox boundaries --- ###
    zgrid = [0, 1500, 60]  # evenly spaced zhalf coords [zmin, zmax, zdelta] [m]
    xgrid = [0, 1500, 50]  # evenly spaced xhalf coords [m]
    ygrid = np.array([0, 100, 200, 300])  # array of yhalf coords [m]

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
    cgrid.write_gridboxboundaries_binary(gridfile, zgrid, xgrid, ygrid, constsfile)
    rgrid.print_domain_info(constsfile, gridfile)

    ### ----- write thermodynamics binaries ----- ###
    thermodyngen = gen_input_thermo.TimeVarying3DThermo(
        PRESSz0, TEMPz0, qvapz0, qcondz0, WMAX, Zlength, Xlength, VMAX, Ylength
    )
    cthermo.write_thermodynamics_binary(
        thermofile, thermodyngen, configfile, constsfile, gridfile
    )

    ### ----- write initial superdroplets binary ----- ###
    nsupers = crdgens.nsupers_at_domain_base(gridfile, constsfile, npergbx, zlim)
    radiigen = rgens.MonoAttrGen(monor)  # all SDs have the same radius [m]
    dryradiigen = dryrgens.ScaledRadiiGen(dryr_sf)  # dryradii are 1/sf of radii [m]
    coord3gen = crdgens.SampleCoordGen(randcoord)  # (not) random coord3 of SDs
    coord1gen = crdgens.SampleCoordGen(randcoord)  # (not) random coord1 of SDs
    coord2gen = crdgens.SampleCoordGen(randcoord)  # (not) random coord2 of SDs
    xiprobdist = probdists.DiracDelta(monor)  # monodisperse droplet probability distrib

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
            constsfile, configfile, gridfile, thermofile, savefigpath, isfigures[1]
        )
        plt.close()
    ### ---------------------------------------------------------------- ###
    ### ---------------------------------------------------------------- ###


if __name__ == "__main__":
    ### args = path2CLEO, path2build, configfile, binpath, gridfile, initSDsfile, thermofile
    main(*sys.argv[1:])
