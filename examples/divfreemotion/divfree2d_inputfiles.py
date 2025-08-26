"""
----- CLEO -----
File: divfree2d_inputfiles.py
Project: divfreemotion
Created Date: Friday 17th November 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
Copyright (c) 2023 MPI-M, Clara Bayley
-----
File Description:
Script generates input files for example of divergence free motion of superdroplets
in a 2-D divergence free wind field.
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
    ### --- essential paths and filenames --- ###
    # path and filenames for creating initial SD conditions
    constants_filename = path2CLEO / "libs" / "cleoconstants.hpp"

    ### --- plotting initialisation figures --- ###
    isfigures = [
        show_figures,
        save_figures,
    ]  # booleans for [showing, saving] initialisation figures
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

    ### --------------------- BINARY FILES GENERATION ---------------------- ###
    ### ----- write gridbox boundaries binary ----- ###
    grid_filename = Path(config["inputfiles"]["grid_filename"])
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
    initsupers_filename = Path(config["initsupers"]["initsupers_filename"])
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
