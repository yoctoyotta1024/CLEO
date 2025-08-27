"""
Copyright (c) 2025 MPI-M, Clara Bayley


----- CLEO -----
File: as2017_inputfiles.py
Project: adiabaticparcel
Created Date: Monday 25th August 2025
Author: Clara Bayley (CB)
Additional Contributors:
-----
Last Modified: Monday 25th August 2025
Modified By: CB
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
Script generates input filesto run adiabatic parcel example which produces plots
similar to Figure 5 of "On the CCN (de)activation nonlinearities" S. Arabas and
S. Shima 2017 to show example of adaibatic parcel expansion and contraction.
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
        "--icond",
        type=int,
        choices=[0, 1, 2],
        default=None,
        help="which initial conditions to generate",
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
    icond=None,
    savefigpath=None,
    show_figures=False,
    save_figures=False,
):
    import numpy as np
    from pathlib import Path
    from ruamel.yaml import YAML

    from cleopy import geninitconds
    from cleopy.initsuperdropsbinary_src import rgens, dryrgens, probdists, attrsgen

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

    # settings for 0D Model (number of SD and grid coordinates)
    nsupers = {0: 64}
    zgrid = np.asarray([0, 100])
    xgrid = np.asarray([0, 100])
    ygrid = np.asarray([0, 100])

    # init_conds dictionary is:
    # icond : [[m^-3] total no. concentration of droplets, monodisperse droplet radius [m]]
    init_conds = {
        0: [500e6, 0.05e-6],
        1: [500e6, 0.1e-6],
        2: [50e6, 0.1e-6],
    }

    # do not generate superdroplet coords
    coord3gen = None
    coord1gen = None
    coord2gen = None

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

    ### ----- write initial superdroplets binary ----- ###
    if gen_supers:
        initsupers_filename = Path(config["initsupers"]["initsupers_filename"])
        numconc, monor = init_conds[icond]

        radiigen = rgens.MonoAttrGen(monor)
        dryradiigen = dryrgens.ScaledRadiiGen(1.0)
        xiprobdist = probdists.DiracDelta(monor)

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
            savelabel=f"_icond{icond}",
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
        icond=args.icond,
        savefigpath=args.savefigpath,
        show_figures=args.show_figures,
        save_figures=args.save_figures,
    )
