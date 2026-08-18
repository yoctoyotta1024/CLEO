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
        "--gen_gbxs",
        action="store_true",  # default is False
        help="Generate gridbox boundaries binary file conditions",
    )
    parser.add_argument(
        "--gen_supers",
        action="store_true",  # default is False
        help="Generate initial superdroplet conditions binary file",
    )
    parser.add_argument(
        "--gen_thermo",
        action="store_true",  # default is False
        help="Generate thermodynamics binary files",
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
    gen_gbxs=False,
    gen_supers=False,
    gen_thermo=False,
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
    pyconfig = config["python_inputfiles"]

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
    # evenly spaced x and z spatial coords [min, max, delta] [m]
    zgrid_spacing = float(pyconfig["zgrid_max"] / pyconfig["zgrid_ngbxs"])
    xgrid_spacing = float(pyconfig["xgrid_max"] / pyconfig["xgrid_ngbxs"])
    zgrid = [0, pyconfig["zgrid_max"], zgrid_spacing]  # evenly spaced zhalf coords [m]
    xgrid = [0, pyconfig["xgrid_max"], xgrid_spacing]  # evenly spaced xhalf coords [m]
    ygrid = np.array([0, pyconfig["ygrid_max"]])  # array of yhalf coords [y0, y1] [m]

    ### --- settings for initial superdroplets --- ###
    # settings for initial superdroplet coordinates
    zlim = pyconfig["sd_zlim"]
    npergbx = pyconfig["nsupers_pergbx"]

    # [min, max] range of initial superdroplet radii (and implicitly solute masses)
    rspan = list(pyconfig["rspan"])

    # settings for initial superdroplet multiplicies (from bimodal Lognormal distribution)
    geomeans = list(pyconfig["geomeans"])
    geosigs = list(pyconfig["geosigs"])
    scalefacs = list(pyconfig["scalefacs"])
    numconc = pyconfig["numconc"]

    ### --- settings for 2D Thermodynamics --- ###
    PRESS = pyconfig["thermo_PRESS"]
    THETA = pyconfig["thermo_THETA"]
    WMAX = pyconfig["thermo_WMAX"]
    Zlength = pyconfig["thermo_Zlength"]
    Xlength = pyconfig["thermo_Xlength"]
    qvapmethod = pyconfig["thermo_qvapmethod"]
    Zbase = pyconfig["thermo_Zbase"]
    sratios = list(pyconfig["thermo_sratios"])
    VVEL = None  # [m/s]
    qcond = 0.0  # background liquid water mass mixing ratio [Kg/Kg]

    ### --------------------- BINARY FILES GENERATION ---------------------- ###
    ### ----- write gridbox boundaries binary ----- ###
    grid_filename = Path(config["inputfiles"]["grid_filename"])
    if gen_gbxs:
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
    if gen_thermo:
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
    if gen_supers:
        initsupers_filename = Path(config["initsupers"]["initsupers_filename"])
        nsupers = crdgens.nsupers_at_domain_base(
            grid_filename, constants_filename, npergbx, zlim
        )
        coord3gen = crdgens.SampleCoordGen(True)  # sample coord3 randomly
        coord1gen = crdgens.SampleCoordGen(True)  # sample coord1 randomly
        coord2gen = None  # do not generate superdroplet coord2s
        xiprobdist = probdists.LnNormal(geomeans, geosigs, scalefacs)
        radiigen = rgens.SampleLog10RadiiGen(
            rspan
        )  # randomly sample radii from rspan [m]
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
        gen_gbxs=args.gen_gbxs,
        gen_supers=args.gen_supers,
        gen_thermo=args.gen_thermo,
        savefigpath=args.savefigpath,
        show_figures=args.show_figures,
        save_figures=args.save_figures,
    )
