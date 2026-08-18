"""
Copyright (c) 2025 MPI-M, Clara Bayley


----- CLEO -----
File: python_bindings_inputfiles.py
Project: python_bindings
Created Date: Friday 6th June 2025
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
Script generates input files for example of using python bindings
(for movement of superdroplets in a 2-D divergence free wind field).
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
    gen_gbxs=False,
    gen_supers=False,
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
    SDgbxs2plt = [
        0
    ]  # gbxindex of initial SDs to plot if any(isfigures) (nb. "all" can be very slow)

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
        gen_gbxs=args.gen_gbxs,
        gen_supers=args.gen_supers,
        savefigpath=args.savefigpath,
        show_figures=args.show_figures,
        save_figures=args.save_figures,
    )
