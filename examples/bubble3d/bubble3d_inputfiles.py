"""
Copyright (c) 2024 MPI-M, Clara Bayley


----- CLEO -----
File: bubble3d_inputfiles.py
Project: bubble3d
Created Date: Friday 17th November 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
Script generates input files for 3D example with time varying
thermodynamics read from ICON output of bubble test case by YAC
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


def get_zgrid(icon_grid_file, num_vertical_levels):
    """returns zgrid for CLEO gridfile with same vertical levels as ICON grid file"""
    import numpy as np
    import xarray as xr

    grid = xr.open_dataset(icon_grid_file)
    idx2 = int(grid.height.values[-1])
    idx1 = int(idx2 - num_vertical_levels - 1)
    zhalf = grid.zghalf.values[idx1:idx2, 0]  # [m]
    zgrid = np.flip(zhalf)

    return zgrid  # [m]


def copy_icon_files(path2build, orginal_icon_grid_file, orginal_icon_data_file):
    import shutil

    assert (
        path2build / "share"
    ).is_dir(), "share directory doesn't exist in build directory"

    icon_grid_file = path2build / "share" / orginal_icon_grid_file.name
    shutil.copyfile(orginal_icon_grid_file, icon_grid_file)

    icon_data_file = path2build / "share" / orginal_icon_data_file.name
    shutil.copyfile(orginal_icon_data_file, icon_data_file)

    return icon_grid_file, icon_data_file


# %%
### -------------------------------- MAIN ---------------------------------- ###
def main(
    path2CLEO,
    path2build,
    config_filename,
    savefigpath=None,
    show_figures=False,
    save_figures=False,
):
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

    ### --- (optional) copy ICON files into build directory for safe-keeping --- ###
    icon_yac_config = config["icon_yac_config"]
    orginal_icon_grid_file = Path(icon_yac_config["orginal_icon_grid_file"])
    orginal_icon_data_file = Path(icon_yac_config["orginal_icon_data_file"])
    icon_grid_file, icon_data_file = copy_icon_files(
        Path(path2build), orginal_icon_grid_file, orginal_icon_data_file
    )

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

    ### --- settings for 3-D gridbox boundaries --- ###
    num_vertical_levels = icon_yac_config["num_vertical_levels"]
    zgrid = get_zgrid(icon_grid_file, num_vertical_levels)  # [m]
    xgrid = [
        0,
        30000,
        2500,
    ]  # evenly spaced xhalf coords [m] # distance must match longitude in config file
    ygrid = [
        0,
        6250,
        1250,
    ]  # evenly spaced xhalf coords [m] # distance must match latitudes in config file

    ### --- settings for initial superdroplets --- ###
    # settings for initial coordinates
    zlim = 1000  # max z coord of superdroplets
    npergbx = 2  # number of superdroplets per gridbox

    # settings for initial radius and aerosol distributions
    monor = 1e-6  # all SDs have this same radius [m]
    dryr_sf = 1.0  # scale factor for dry radii [m]
    numconc = 5e8  # total no. conc of real droplets [m^-3]
    randcoord = False  # sample SD spatial coordinates randomly or not

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

    ### ----- write initial superdroplets binary ----- ###
    initsupers_filename = Path(config["initsupers"]["initsupers_filename"])
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
        savefigpath=args.savefigpath,
        show_figures=args.show_figures,
        save_figures=args.save_figures,
    )
