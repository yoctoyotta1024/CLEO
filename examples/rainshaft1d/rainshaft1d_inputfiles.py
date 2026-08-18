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
    pyconfig = config["python_inputfiles"]

    ### ------------------------ INPUT PARAMETERS -------------------------- ###
    ### --- required CLEO cleoconstants.hpp file --- ###
    constants_filename = Path(config["inputfiles"]["constants_filename"])

    ### --- plots of initial conditions --- ###
    isfigures = [
        show_figures,
        save_figures,
    ]  # booleans for [showing, saving] initialisation figures
    SDgbxs2plt = list(
        range(min(39, config["domain"]["ngbxs"]), config["domain"]["ngbxs"])
    )  # gbxindex of initial SDs to plot if any(isfigures) (nb. "all" can be very slow)
    SDgbxs2plt = [random.choice(SDgbxs2plt)]  # choose random gbx from list to plot

    ### --- settings for 1-D gridbox boundaries --- ###
    # evenly spaced z spatial coords [min, max, delta] [m]
    zgrid_spacing = float(pyconfig["zgrid_max"] / config["domain"]["ngbxs"])
    zgrid = [0, pyconfig["zgrid_max"], zgrid_spacing]  # evenly spaced zhalf coords [m]
    xgrid = np.array([0, pyconfig["xgrid_max"]])  # array of xhalf coords [x0, x1] [m]
    ygrid = np.array([0, pyconfig["ygrid_max"]])  # array of yhalf coords [y0, y1] [m]

    ### --- settings for 1-D Thermodynamics --- ###
    PRESSz0 = pyconfig["thermo_PRESSz0"]  # [Pa]
    TEMPz0 = pyconfig["thermo_TEMPz0"]
    qvapz0 = pyconfig["thermo_qvapz0"]
    Zbase = pyconfig["thermo_Zbase"]
    TEMPlapses = list(pyconfig["thermo_TEMPlapses"])
    qvaplapses = list(pyconfig["thermo_qvaplapses"])
    WVEL = pyconfig["thermo_WVEL"]
    Wlength = pyconfig["thermo_Wlength"]
    qcond = 0.0  # background liquid water mass mixing ratio [Kg/Kg]

    ### --- settings for initial superdroplets --- ###
    # initial superdroplet coordinates
    zlim = pyconfig["sd_zlim"]
    npergbx = pyconfig["nsupers_pergbx"]

    # initial superdroplet radii (and implicitly solute masses)
    rspan = list(pyconfig["rspan"])
    dryr_sf = pyconfig["dryr_sf"]

    # settings for initial superdroplet multiplicies (from trimodal Lognormal distribution)
    geomeans = list(pyconfig["geomeans"])
    geosigs = list(pyconfig["geosigs"])
    scalefacs = list(pyconfig["scalefacs"])
    numconc = pyconfig["numconc"]

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
            isprintinfo=True,
            isfigures=isfigures,
            savefigpath=savefigpath,
        )

    ### ----- write thermodynamics binaries ----- ###
    if gen_thermo:
        thermog = thermogen.HydrostaticLapseRates(
            config_filename,
            constants_filename,
            PRESSz0,
            TEMPz0,
            qvapz0,
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
    if gen_supers:
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
        gen_gbxs=args.gen_gbxs,
        gen_supers=args.gen_supers,
        gen_thermo=args.gen_thermo,
        savefigpath=args.savefigpath,
        show_figures=args.show_figures,
        save_figures=args.save_figures,
    )
