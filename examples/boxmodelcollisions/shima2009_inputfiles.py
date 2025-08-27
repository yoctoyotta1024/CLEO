"""
Copyright (c) 2025 MPI-M, Clara Bayley


----- CLEO -----
File: shima2009_inputfiles.py
Project: boxmodelcollisions
Created Date: Thursday 21st August 2025
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
Script for generating input files for CLEO 0-D box model for collisions using the
Golovin or Long (hydrodynamic) kernel in a way comparable to Shima et al. 2009 Fig. 2
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
        "kernel",
        type=str,
        choices=["golovin", "long1", "long2"],
        help="kernel example to run",
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
### ------------------------- FUNCTION DEFINITIONS ------------------------- ###
def initial_superdroplet_conditions_for_setup(
    path2CLEO,
    config_filename,
    nsupers,
    radiigen,
    xiprobdist,
    numconc,
    isfigures=[False, False],
    savefigpath=None,
    savelabel="",
    gbxs2plt="all",
):
    from pathlib import Path
    from ruamel.yaml import YAML

    from cleopy import geninitconds
    from cleopy.initsuperdropsbinary_src import rgens, attrsgen

    ### --- Load the config YAML file --- ###
    yaml = YAML()
    with open(config_filename, "r") as file:
        config = yaml.load(file)

    ### --- settings for superdroplet attributes --- ###
    dryradius = 1e-16  # all SDs have negligible solute [m]
    dryradiigen = rgens.MonoAttrGen(dryradius)
    coord3gen = None  # do not generate superdroplet coords
    coord1gen = None
    coord2gen = None

    initattrsgen = attrsgen.AttrsGenerator(
        radiigen, dryradiigen, xiprobdist, coord3gen, coord1gen, coord2gen
    )

    ### ----- write initial superdroplets binary ----- ###
    constants_filename = Path(config["inputfiles"]["constants_filename"])
    initsupers_filename = Path(config["initsupers"]["initsupers_filename"])
    grid_filename = Path(config["inputfiles"]["grid_filename"])
    geninitconds.generate_initial_superdroplet_conditions(
        initattrsgen,
        initsupers_filename,
        config_filename,
        constants_filename,
        grid_filename,
        nsupers,
        numconc,
        isprintinfo=True,
        isfigures=isfigures,
        savefigpath=savefigpath,
        gbxs2plt=gbxs2plt,
        savelabel=savelabel,
    )


# %%
### -------------------------------- MAIN ---------------------------------- ###
def main(
    path2CLEO,
    path2build,
    config_filename,
    kernel,
    savefigpath=None,
    show_figures=False,
    save_figures=False,
):
    import numpy as np
    from pathlib import Path
    from ruamel.yaml import YAML

    import attrgens_shima2009
    from cleopy import geninitconds

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
    SDgbxs2plt = "all"  # gbxindex of initial SDs to plot if any(isfigures) (nb. "all" can be very slow)

    ### --- settings for 0-D Model gridbox boundaries --- ###
    zgrid = np.asarray([0, 100])
    xgrid = np.asarray([0, 100])
    ygrid = np.asarray([0, 100])

    ### --- settings for initial superdroplets --- ###
    nsupers = int(config["domain"]["maxnsupers"])

    ### --- settings for initial superdroplets for setup 1 --- ###
    # radius distirbution from exponential in droplet volume
    rspan_1 = [0.62e-6, 6.34e-2]  # max and min range of radii to sample [m]
    volexpr0_1 = 30.531e-6  # peak of volume exponential distribution [m]
    numconc_1 = 2**23  # total no. conc of real droplets [m^-3]

    # attribute generators
    xiprobdist_1 = attrgens_shima2009.SampleXiShima2009()
    radiigen_1 = attrgens_shima2009.SampleRadiiShima2009(
        volexpr0_1, rspan_1
    )  # radii are sampled from rspan [m]

    ### --- settings for initial superdroplets for setup 2 --- ###
    # radius distirbution from exponential in droplet volume
    rspan_2 = [0.62e-6, 6.34e-2]  # max and min range of radii to sample [m]
    volexpr0_2 = 10.117e-6  # peak of volume exponential distribution [m]
    numconc_2 = 27 * 2**23  # total no. conc of real droplets [m^-3]

    # attribute generators
    xiprobdist_2 = attrgens_shima2009.SampleXiShima2009()
    radiigen_2 = attrgens_shima2009.SampleRadiiShima2009(
        volexpr0_2, rspan_2
    )  # radii are sampled from rspan [m]

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

    ### ----- write initial superdroplets binary ----- ###
    if "golovin" == kernel or "long1" == kernel:
        initial_superdroplet_conditions_for_setup(
            path2CLEO,
            config_filename,
            nsupers,
            radiigen_1,
            xiprobdist_1,
            numconc_1,
            isfigures=isfigures,
            savefigpath=savefigpath,
            savelabel=f"_{kernel}",
            gbxs2plt=SDgbxs2plt,
        )
    elif "long2" == kernel:
        initial_superdroplet_conditions_for_setup(
            path2CLEO,
            config_filename,
            nsupers,
            radiigen_2,
            xiprobdist_2,
            numconc_2,
            isfigures=isfigures,
            savefigpath=savefigpath,
            savelabel=f"_{kernel}",
            gbxs2plt=SDgbxs2plt,
        )
    else:
        raise ValueError(
            "kernel for examples not recognised, please choose from: golovin, long1 or long2"
        )


# %%
### --------------------------- RUN PROGRAM -------------------------------- ###
if __name__ == "__main__":
    args = parse_arguments()
    main(
        args.path2CLEO,
        args.path2build,
        args.config_filename,
        args.kernel,
        savefigpath=args.savefigpath,
        show_figures=args.show_figures,
        save_figures=args.save_figures,
    )
