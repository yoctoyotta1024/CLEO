"""
Copyright (c) 2025 MPI-M, Clara Bayley


----- CLEO -----
File: breakup_inputfiles.py
Project: boxmodelcollisions
Created Date: Friday 22nd August 2025
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
Script for generating input files for CLEO 0-D box model for collisions using
selected collision kernels with breakup (e.g. Low and Lists's) in a way
comparable to Shima et al. 2009 Fig. 2(b)
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
        choices=["long", "lowlist", "szakallurbich", "testikstraub"],
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

    from cleopy import geninitconds
    from cleopy.initsuperdropsbinary_src import rgens, probdists, attrsgen

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
    # settings for superdroplet coordinates
    nsupers = int(config["domain"]["maxnsupers"])

    # settings for superdroplet attributes
    dryradius = 1e-16  # all SDs have negligible solute [m]
    coord3gen = None  # do not generate superdroplet coords
    coord1gen = None
    coord2gen = None

    # radius distirbution from exponential in droplet volume for setup 1
    rspan = [1e-7, 9e-5]  # max and min range of radii to sample [m]
    volexpr0 = 30.531e-6  # peak of volume exponential distribution [m]
    numconc = 2**23  # total no. conc of real droplets [m^-3]

    # attribute generators
    xiprobdist = probdists.VolExponential(volexpr0, rspan)
    radiigen = rgens.SampleLog10RadiiGen(rspan)  # radii are sampled from rspan [m]
    dryradiigen = rgens.MonoAttrGen(dryradius)

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
    initsupers_filename = Path(config["initsupers"]["initsupers_filename"])
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
        isprintinfo=True,
        isfigures=isfigures,
        savefigpath=savefigpath,
        gbxs2plt=SDgbxs2plt,
        savelabel=f"_{kernel}",
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
