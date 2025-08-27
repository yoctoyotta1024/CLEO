"""
Copyright (c) 2025 MPI-M, Clara Bayley


----- CLEO -----
File: rainshaft1d_inputfiles.py
Project: rainshaft1d
Created Date: Friday 22nd August 2025
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
Script generates input files, for example of 1-D rainshaft with constant
thermodynamics read from a file.
"""


# %%
### ------------------------- FUNCTION DEFINITIONS ------------------------- ###
def parse_arguments():
    import argparse
    from pathlib import Path

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "path2CLEO", type=Path, help="Absolute path to CLEO directory (for cleopy)"
    )
    parser.add_argument(
        "path2build", type=Path, help="Absolute path to build directory"
    )
    parser.add_argument(
        "config_filename", type=Path, help="Absolute path to configuration YAML file"
    )
    parser.add_argument(
        "thermofiles",
        type=Path,
        help="Absolute path to derive thermoynamics binary files",
    )
    parser.add_argument(
        "--savefigpath",
        type=Path,
        default=None,
        help="Directory to save initialiation figures in (is save_figures is True)",
    )
    parser.add_argument(
        "--show_figures",
        action="store_true",  # default is False
        help="Show initialiation figures",
    )
    parser.add_argument(
        "--save_figures",
        action="store_true",  # default is False
        help="Save initialiation figures in savefigpath",
    )
    return parser.parse_args()


# %%
### -------------------------------- MAIN ---------------------------------- ###
def main(
    path2CLEO,
    path2build,
    config_filename,
    thermofiles,
    savefigpath=None,
    show_figures=False,
    save_figures=False,
):
    import numpy as np
    import random
    from pathlib import Path
    from ruamel.yaml import YAML

    from cleopy import geninitconds
    from cleopy.initsuperdropsbinary_src import (
        crdgens,
        rgens,
        dryrgens,
        probdists,
        attrsgen,
    )
    from cleopy.thermobinary_src import thermogen, windsgen, thermodyngen

    if path2CLEO == path2build:
        raise ValueError("build directory cannot be CLEO")

    ### --- Load the config YAML file --- ###
    yaml = YAML()
    with open(config_filename, "r") as file:
        config = yaml.load(file)

    ### ------------------------ INPUT PARAMETERS -------------------------- ###
    ### --- required CLEO cleoconstants.hpp file --- ###
    constants_filename = Path(config["inputfiles"]["constants_filename"])

    ### --- plots of initial conditions --- ###
    isfigures = [
        show_figures,
        save_figures,
    ]  # booleans for [showing, saving] initialisation figures
    SDgbxs2plt = list(
        range(39, 124)
    )  # gbxindex of initial SDs to plot if any(isfigures) (nb. "all" can be very slow)
    SDgbxs2plt = [random.choice(SDgbxs2plt)]  # choose random gbx from list to plot

    ### --- settings for 1-D gridbox boundaries --- ###
    zgrid = [0, 2500, 20]  # evenly spaced zhalf coords [zmin, zmax, zdelta] [m]
    xgrid = np.array([0, 20])  # array of xhalf coords [m]
    ygrid = np.array([0, 20])  # array of yhalf coords [m]

    ### --- settings for 1-D Thermodynamics --- ###
    PRESS0 = 101315  # [Pa]
    TEMP0 = 297.9  # [K]
    qvap0 = 0.016  # [Kg/Kg]
    Zbase = 800  # [m]
    TEMPlapses = [9.8, 6.5]  # -dT/dz [K/km]
    qvaplapses = [2.97, "saturated"]  # -dvap/dz [g/Kg km^-1]
    qcond = 0.0  # [Kg/Kg]
    WVEL = 4.0  # [m/s]
    Wlength = 1000  # [m] use constant W (Wlength=0.0), or sinusoidal 1-D profile below cloud base

    ### --- settings for initial superdroplets --- ###
    # initial superdroplet coordinates
    zlim = 800  # min z coord of superdroplets [m]
    npergbx = 256  # number of superdroplets per gridbox

    # initial superdroplet radii (and implicitly solute masses)
    rspan = [3e-9, 5e-5]  # min and max range of radii to sample [m]
    dryr_sf = 1.0  # dryradii are 1/sf of radii [m]

    # settings for initial superdroplet multiplicies
    geomeans = [0.02e-6, 0.2e-6, 3.5e-6]
    geosigs = [1.55, 2.3, 2]
    scalefacs = [1e6, 0.3e6, 0.025e6]
    numconc = np.sum(scalefacs) * 1000

    ### --------------------- BINARY FILES GENERATION ---------------------- ###
    ### ----- write gridbox boundaries binary ----- ###
    grid_filename = Path(config["inputfiles"]["grid_filename"])
    geninitconds.generate_gridbox_boundaries(
        grid_filename,
        zgrid,
        xgrid,
        ygrid,
        constants_filename,
        isprintinfo=True,
        isfigures=isfigures,
        savefigpath=savefigpath,
    )

    ### ----- write thermodynamics binaries ----- ###
    thermog = thermogen.HydrostaticLapseRates(
        config_filename,
        constants_filename,
        PRESS0,
        TEMP0,
        qvap0,
        Zbase,
        TEMPlapses,
        qvaplapses,
        qcond,
    )
    windsg = windsgen.SinusoidalUpdraught(WVEL, None, None, Wlength)
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
    initsupers_filename = Path(config["initsupers"]["initsupers_filename"])
    nsupers = crdgens.nsupers_at_domain_top(
        grid_filename, constants_filename, npergbx, zlim
    )
    coord3gen = crdgens.SampleCoordGen(True)  # sample coord3 randomly
    coord1gen = None  # do not generate superdroplet coord2s
    coord2gen = None  # do not generate superdroplet coord2s

    xiprobdist = probdists.LnNormal(geomeans, geosigs, scalefacs)
    radiigen = rgens.SampleLog10RadiiGen(rspan)
    dryradiigen = dryrgens.ScaledRadiiGen(dryr_sf)

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


# %%
### --------------------------- RUN PROGRAM -------------------------------- ###
if __name__ == "__main__":
    args = parse_arguments()
    main(
        args.path2CLEO,
        args.path2build,
        args.config_filename,
        args.thermofiles,
        savefigpath=args.savefigpath,
        show_figures=args.show_figures,
        save_figures=args.save_figures,
    )
