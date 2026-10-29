"""
Copyright (c) 2026 MPI-M, Clara Bayley


----- CLEO -----
File: roughpaper_inputfiles.py
Project: roughpaper
Created Date: Wednesday 7th October 2026
Author: Clara Bayley (CB)
Additional Contributors: Harshada Balasubramanian (HB)
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
Script generates input files for the roughpaper executable "cleocoupledsdm". It is also an
example of various ways to use the cleopy module to create the binary files for the gridbox
boundaries, the initial superdroplet conditions and the thermodynamics (when using coupled
dynamics "fromfile") to read into CLEO SDM. (Previously the scripts
create_gbxboundariesbinary_script.py, create_initsuperdropsbinary_script.py and
create_thermobinaries_script.py.)
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


def generate_gridbox_boundaries(
    grid_filename, constants_filename, isfigures, savefigpath
):
    """example of using cleopy to create the gridbox boundaries binary file"""
    import numpy as np
    from cleopy import geninitconds

    ### input parameters for zcoords of gridbox boundaries
    zmax = 100  # maximum z coord [m]
    zmin = 0  # minimum z coord [m]
    zdelta = 10  # even spacing
    zgrid = [zmin, zmax, zdelta]
    # zgrid = np.arange(zmin, zmax+zdelta, zdelta)

    ### input parameters for x coords of gridbox boundaries
    xgrid = [0, 20, 20]

    ### input parameters for y coords of gridbox boundaries
    ygrid = np.asarray([0, 20])

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


def generate_thermodynamics(
    thermofiles,
    config_filename,
    constants_filename,
    grid_filename,
    isfigures,
    savefigpath,
):
    """example of using cleopy to create the thermodynamics binary files
    (for coupled dynamics read from file)"""
    from cleopy import geninitconds
    from cleopy.thermobinary_src import thermogen, windsgen, thermodyngen

    ### --- Choose Initial Thermodynamic Conditions for Gridboxes  --- ###

    ### --- Thermo (temp, press, qvap and cond) Conditions  --- ###

    # ### --- Constant and Uniform --- ###
    # P_INIT = 101500.0                       # initial pressure [Pa]
    # TEMP_INIT = 288.15                      # initial parcel temperature [T]
    # relh_init = 0.999                       # initial relative humidity (%)
    # qvap = None                             # use relative humidity to set qvap
    # qc_init = 0.0                           # initial liquid water content []
    # thermog = thermogen.ConstUniformThermo(P_INIT, TEMP_INIT, None,
    #                              qc_init, relh=relh_init,
    #                              constants_filename=constants_filename)

    ### --- 1-D T and qv set by Lapse Rates --- ###
    PRESSz0 = 101315  # [Pa]
    TEMPz0 = 297.9  # [K]
    qvapz0 = 0.016  # [Kg/Kg]
    Zbase = 800  # [m]
    TEMPlapses = [9.8, 6.5]  # -dT/dz [K/km]
    qvaplapses = [2.97, "saturated"]  # -dvap/dz [g/Kg km^-1]
    qcond = 0.0  # [Kg/Kg]
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

    # ### --- Hydrostatic Dry Adiabat --- ###
    # ### ---   or Simple z Profile   --- ###
    # PRESSz0 = 101500  # [Pa]
    # THETA = 289  # [K]
    # qcond = 0.0  # [Kg/Kg]
    # qvapmethod = "sratio"
    # Zbase = 750  # [m]
    # sratios = [0.85, 1.0001]  # s_ratio [below, above] Zbase
    # # moistlayer = False
    # moistlayer = {"z1": 700, "z2": 800, "x1": 0, "x2": 750, "mlsratio": 1.005}
    # thermog = thermogen.DryHydrostaticAdiabatic2TierRelH(
    #     config_filename,
    #     constants_filename,
    #     PRESSz0,
    #     THETA,
    #     qvapmethod,
    #     sratios,
    #     Zbase,
    #     qcond,
    #     moistlayer,
    # )
    # thermog = thermogen.Simple2TierRelativeHumidity(config_filename, constants_filename, PRESSz0,
    #                                         THETA, qvapmethod, sratios, Zbase,
    #                                         qcond)

    ### --- Wind Field (wvel, vvel, uvel) Conditions  --- ###

    ### --- Constant and Uniform --- ###
    W_INIT = 0.0  # initial vertical (coord3) velocity [m/s]
    U_INIT = 0.0  # initial eastwards (coord1) velocity [m/s]
    V_INIT = 0.0  # initial northwards (coord2) velocity [m/s]
    windsg = windsgen.ConstUniformWinds(W_INIT, U_INIT, V_INIT)

    # ### --- 1D Vertical Sinusoid --- ###
    # WMAX = 0.0  # [m/s]
    # UVEL = None
    # VVEL = None
    # Wlength = (
    #     1000  # [m] use constant W (Wlength=0.0), or sinusoidal 1-D profile below cloud base
    # )
    # windsg = windsgen.SinusoidalUpdraught(WMAX, UVEL, VVEL, Wlength)

    # ### --- 2D Flow Field --- ###
    # WMAX = 0.6  # [m/s]
    # VVEL = 1.0  # [m/s]
    # Zlength = 1500  # [m]
    # Xlength = 1500  # [m]
    # # windsg = thermog.create_default_windsgen(
    # #     WMAX, Zlength, Xlength, VVEL
    # # )  # only for DryHydrostaticAdiabatic2TierRelH
    # windsg = windsgen.Simple2DFlowField(
    #     config_filename, constants_filename, WMAX, Zlength, Xlength, VVEL
    # )

    ### --- Thermodynamic + Winds Conditions  --- ###
    thermodyngenerator = thermodyngen.ThermodynamicsGenerator(thermog, windsg)

    geninitconds.generate_thermodynamics_conditions_fromfile(
        thermofiles,
        thermodyngenerator,
        config_filename,
        constants_filename,
        grid_filename,
        isfigures=isfigures,
        savefigpath=savefigpath,
    )


def generate_initial_superdroplets(
    initsupers_filename,
    config_filename,
    constants_filename,
    grid_filename,
    isfigures,
    savefigpath,
):
    """example of various ways to use cleopy to create the binary file for the initial
    superdroplet conditions"""
    # import numpy as np
    from cleopy import geninitconds
    from cleopy.initsuperdropsbinary_src import rgens, probdists, attrsgen, crdgens

    gbxs2plt = "all"  # indexes of GBx index of SDs to plot (nb. "all" can be very slow)

    ### --- Number of Superdroplets per Gridbox --- ###
    ### ---        (an int or dict of ints)     --- ###
    # zlim = 800
    # npergbx = 8192
    # nsupers =  crdgens.nsupers_at_domain_base(grid_filename, constants_filename, npergbx, zlim) # supers where z <= zlim
    # nsupers = crdgens.nsupers_at_domain_top(
    #     grid_filename, constants_filename, npergbx, zlim
    # )  # supers where z >= zlim
    nsupers = 256
    ### ------------------------------------------- ###

    ### --- Choice of Superdroplet Radii Generator --- ###
    # monor                = 0.05e-6                        # all SDs have this same radius [m]
    # radiigen  =  rgens.MonoAttrGen(monor)                 # all SDs have the same radius [m]

    rspan = [1e-6, 1e-3]  # min and max range of radii to sample [m]
    radiigen = rgens.SampleLog10RadiiGen(rspan)  # radii are sampled from rspan [m]
    ### ---------------------------------------------- ###

    ### --- Choice of Superdroplet Dry Radii Generator --- ###
    monodryr = 1e-9  # all SDs have this same dryradius [m]
    dryradiigen = rgens.MonoAttrGen(monodryr)  # all SDs have the same dryradius [m]

    # dryr_sf = 1.0  # scale factor for dry radii [m]
    # dryradiigen = dryrgens.ScaledRadiiGen(dryr_sf)  # dryradii are 1/sf of radii [m]

    ### ---------------------------------------------- ###

    ### --- Choice of Droplet Radius Probability Distribution --- ###
    # dirac0               = monor                         # radius in sample closest to this value is dirac delta peak
    # numconc              = 1e6                         # total no. conc of real droplets [m^-3]
    # numconc              = 512e6                         # total no. conc of real droplets [m^-3]
    # xiprobdist = probdists.DiracDelta(dirac0)

    # geomeans           = [0.075e-6]                  # lnnormal modes' geometric mean droplet radius [m]
    # geosigs            = [1.5]                       # lnnormal modes' geometric standard deviation
    # scalefacs          = [1]                         # relative heights of modes
    # geomeans             = [0.02e-6, 0.2e-6, 3.5e-6]
    # geosigs              = [1.55, 2.3, 2]
    # scalefacs            = [1, 0.3, 0.025]
    # geomeans = [0.02e-6, 0.15e-6]
    # geosigs = [1.4, 1.6]
    # scalefacs = [0.6, 0.4]
    # numconc = np.sum(scalefacs) * 1e9
    # xiprobdist = probdists.LnNormal(geomeans, geosigs, scalefacs)

    # volexpr0             = 30.531e-6                   # peak of volume exponential distribution [m]
    # numconc              = 2**(23)                     # total no. conc of real droplets [m^-3]
    # xiprobdist = probdists.VolExponential(volexpr0, rspan)

    reff = 7e-6  # effective radius [m]
    nueff = 0.08  # effective variance
    # xiprobdist = probdists.ClouddropsHansenGamma(reff, nueff)
    rdist1 = probdists.ClouddropsHansenGamma(reff, nueff)
    nrain = 3000  # raindrop concentration [m^-3]
    qrain = 0.9  # rainwater content [g/m^3]
    dvol = 8e-4  # mean volume diameter [m]
    # xiprobdist = probdists.RaindropsGeoffroyGamma(nrain, qrain, dvol)
    rdist2 = probdists.RaindropsGeoffroyGamma(nrain, qrain, dvol)
    numconc = 1e9  # [m^3]
    distribs = [rdist1, rdist2]
    scalefacs = [1000, 1]
    xiprobdist = probdists.CombinedRadiiProbDistribs(distribs, scalefacs)

    ### --------------------------------------------------------- ###

    ### --- Choice of Superdroplet Coord3 Generator --- ###
    # monocoord3           = 1000                        # all SDs have this same coord3 [m]
    # coord3gen            =  crdgens.MonoCoordGen(monocoord3)
    coord3gen = crdgens.SampleCoordGen(True)  # sample coord3 range randomly or not
    # coord3gen = None  # do not generate superdroplet coord3s
    ### ----------------------------------------------- ###

    ### --- Choice of Superdroplet Coord1 Generator --- ###
    # monocoord1           = 200                        # all SDs have this same coord1 [m]
    # coord1gen            =  crdgens.MonoCoordGen(monocoord1)
    # coord1gen = crdgens.SampleCoordGen(True)  # sample coord1 range randomly or not
    coord1gen = None  # do not generate superdroplet coord1s
    ### ----------------------------------------------- ###

    ### --- Choice of Superdroplet Coord2 Generator --- ###
    # monocoord2           = 1000                        # all SDs have this same coord2 [m]
    # coord2gen            =  crdgens.MonoCoordGen(monocoord2)
    # coord2gen = crdgens.SampleCoordGen(True)  # sample coord1 range randomly or not
    coord2gen = None  # do not generate superdroplet coord2s
    ### ----------------------------------------------- ###

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
        gbxs2plt=gbxs2plt,
    )


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
    from pathlib import Path
    from ruamel.yaml import YAML

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

    ### --------------------- BINARY FILES GENERATION ---------------------- ###
    ### ----- write gridbox boundaries binary ----- ###
    grid_filename = Path(config["inputfiles"]["grid_filename"])
    if gen_gbxs:
        generate_gridbox_boundaries(
            grid_filename, constants_filename, isfigures, savefigpath
        )

    ### ----- write thermodynamics binaries ----- ###
    if gen_thermo:
        generate_thermodynamics(
            thermofiles,
            config_filename,
            constants_filename,
            grid_filename,
            isfigures,
            savefigpath,
        )

    ### ----- write initial superdroplets binary ----- ###
    if gen_supers:
        initsupers_filename = Path(config["initsupers"]["initsupers_filename"])
        generate_initial_superdroplets(
            initsupers_filename,
            config_filename,
            constants_filename,
            grid_filename,
            isfigures,
            savefigpath,
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
