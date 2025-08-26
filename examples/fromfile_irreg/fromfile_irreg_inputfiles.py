"""
Copyright (c) 2024 MPI-M, Clara Bayley


----- CLEO -----
File: fromfile_irreg_inputfiles.py
Project: fromfile_irreg
Created Date: Wednesday 11th September 2024
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
Script generates input files for 3D example with irregular grid and
time varying thermodynamics read from binary files.
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

    from src import gen_input_thermo
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

    ### ------------------------ INPUT PARAMETERS -------------------------- ###
    ### --- required CLEO cleoconstants.hpp file --- ###
    constants_filename = Path(config["inputfiles"]["constants_filename"])

    ### --- booleans for [showing, saving] initialisation figures --- ###
    isfigures = [show_figures, save_figures]

    ### --- settings for 3-D irregular gridbox boundaries --- ###
    zgrid = np.array([0, 20, 30, 45, 60, 80, 90, 120, 140, 180, 360, 500, 1000, 1500])
    xgrid = np.array(
        [0, 33, 205, 440, 650, 915, 1033, 1100, 1300, 1450, 1500]
    )  # evenly spaced xhalf coords [m]
    ygrid = np.array([0, 10, 75, 100, 150, 200, 300])  # array of yhalf coords [m]

    ### --- settings for initial superdroplets --- ###
    # settings for initial superdroplet coordinates
    zlim = 1000  # max z coord of superdroplets
    npergbx = 2  # number of superdroplets per gridbox

    # settings for initial radius and aerosol distributions
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
    thermodyngen = gen_input_thermo.TimeVarying3DThermodyn(
        PRESSz0, TEMPz0, qvapz0, qcondz0, WMAX, Zlength, Xlength, VMAX, Ylength
    )
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
