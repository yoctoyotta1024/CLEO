"""
Copyright (c) 2024 MPI-M, Clara Bayley


----- CLEO -----
File: geninitconds.py
Project: cleopy
Created Date: Friday 8th November 2024
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
Script to help with generating and plotting input files
"""


import matplotlib.pyplot as plt


def generate_gridbox_boundaries(
    grid_filename,
    zgrid,
    xgrid,
    ygrid,
    constants_filename,
    isprintinfo=False,
    isfigures=[False, False],
    savefigpath=None,
    savelabel="",
):
    """write gridbox boundaries binary. [shows, saves] figures if isfigures = [True, True]"""
    from .gbxboundariesbinary_src import create_gbxboundaries as cgrid
    from .gbxboundariesbinary_src import read_gbxboundaries as rgrid

    cgrid.write_gridboxboundaries_binary(
        grid_filename, zgrid, xgrid, ygrid, constants_filename
    )

    if isprintinfo:
        rgrid.print_domain_info(constants_filename, grid_filename)

    if any(isfigures):
        rgrid.plot_gridboxboundaries(
            constants_filename,
            grid_filename,
            savefig=isfigures[1],
            savefigpath=savefigpath,
            savelabel=savelabel,
        )
        if isfigures[0]:
            plt.show()
        plt.close("all")


def generate_initial_superdroplet_conditions(
    initattrsgen,
    initsupers_filename,
    config_filename,
    constants_filename,
    grid_filename,
    nsupers,
    numconc,
    numconc_tolerance=0.0,
    isprintinfo=False,
    isfigures=[False, False],
    savefigpath=None,
    gbxs2plt=0,
    savelabel="",
):
    """write initial superdroplets binary. [shows, saves] figures if isfigures = [True, True]"""
    from .initsuperdropsbinary_src import create_initsuperdrops as csupers
    from .initsuperdropsbinary_src import read_initsuperdrops as rsupers

    csupers.write_initsuperdrops_binary(
        initsupers_filename,
        initattrsgen,
        config_filename,
        constants_filename,
        grid_filename,
        nsupers,
        numconc,
        numconc_tolerance=numconc_tolerance,
        isprint=isprintinfo,
    )

    if isprintinfo:
        rsupers.print_initsupers_infos(
            initsupers_filename, config_filename, constants_filename, grid_filename
        )

    if any(isfigures):
        rsupers.plot_initGBxs_distribs(
            config_filename,
            constants_filename,
            initsupers_filename,
            grid_filename,
            gbxs2plt,
            savefig=isfigures[1],
            savefigpath=savefigpath,
            savelabel=savelabel,
        )
        if isfigures[0]:
            plt.show()
        plt.close("all")


def generate_thermodynamics_conditions_fromfile(
    thermofiles,
    thermodyngen,
    config_filename,
    constants_filename,
    grid_filename,
    isfigures=[False, False],
    savefigpath=None,
    savelabel="",
):
    """write thermodynamics binaries. [shows, saves] figures if isfigures = [True, True]"""
    from .thermobinary_src import create_thermodynamics as cthermo
    from .thermobinary_src import read_thermodynamics as rthermo

    cthermo.write_thermodynamics_binary(
        thermofiles,
        thermodyngen,
        config_filename,
        constants_filename,
        grid_filename,
    )

    if any(isfigures):
        rthermo.plot_thermodynamics(
            constants_filename,
            config_filename,
            grid_filename,
            thermofiles,
            savefig=isfigures[1],
            savefigpath=savefigpath,
            savelabel=savelabel,
        )
        if isfigures[0]:
            plt.show()
        plt.close("all")
