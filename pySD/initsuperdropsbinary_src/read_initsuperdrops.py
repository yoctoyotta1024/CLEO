"""
Copyright (c) 2024 MPI-M, Clara Bayley


----- CLEO -----
File: read_initsuperdrops.py
Project: initsuperdropsbinary_src
Created Date: Friday 13th October 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
"""

import numpy as np
import matplotlib.pyplot as plt

from .create_initsuperdrops import initsupers_inputsdict, ManyAttrs
from ..readbinary import readbinary
from ..gbxboundariesbinary_src.read_gbxboundaries import get_gbxvols_from_gridfile


def plot_initGBxs_distribs(
    config_filename,
    constants_filename,
    initsupers_filename,
    grid_filename,
    gbxs2plt,
    savefig=False,
    savefigpath=None,
    savelabel="",
):
    """plot initial superdroplet distribution from initsupersfile binary
    of every gridbox with index in gbx2plts"""

    plot_initGBxs_attrdistribs(
        config_filename,
        constants_filename,
        initsupers_filename,
        grid_filename,
        gbxs2plt,
        savefigpath,
        savefig,
        savelabel,
    )
    plot_initGBxs_dropletmasses(
        config_filename,
        constants_filename,
        initsupers_filename,
        grid_filename,
        gbxs2plt,
        savefigpath,
        savefig,
        savelabel,
    )


def get_superdroplet_attributes(
    config_filename, constants_filename, initsupers_filename
):
    """get gridbox boundaries from binary file and
    re-dimensionalise usign COORD0 const from constants_filename"""

    inputs = initsupers_inputsdict(config_filename, constants_filename)

    attrs = read_dimless_superdrops_binary(initsupers_filename, isprint=False)

    # re-dimensionalise SD attributes
    attrs.radius = attrs.radius * inputs["R0"]
    attrs.msol = attrs.msol * inputs["MASS0"]
    attrs.coord3 = attrs.coord3 * inputs["COORD0"]
    attrs.coord1 = attrs.coord1 * inputs["COORD0"]
    attrs.coord2 = attrs.coord2 * inputs["COORD0"]

    return attrs


def read_dimless_superdrops_binary(filename, isprint=True):
    """return dimenionsless gbx boundaries by reading binary file"""

    datatypes = [np.uintc, np.uint, np.double, np.double]
    datatypes += [np.double] * 3
    data, ndata_pervar = readbinary(filename, isprint=isprint)

    ll = [0, 0, 0, 0, 0, 0]  # indexs for division of data list between each variable
    for n in range(1, len(ndata_pervar)):
        ll[n - 1] = np.sum(ndata_pervar[:n])

    attrs = ManyAttrs()
    attrs.sdgbxindex = np.asarray(data[: ll[0]], dtype=datatypes[0])
    attrs.xi = np.asarray(data[ll[0] : ll[1]], dtype=datatypes[1])
    attrs.radius = np.asarray(data[ll[1] : ll[2]], dtype=datatypes[2])
    attrs.msol = np.asarray(data[ll[2] : ll[3]], dtype=datatypes[3])
    attrs.coord3 = np.asarray(data[ll[3] : ll[4]], dtype=datatypes[4])
    attrs.coord1 = np.asarray(data[ll[4] : ll[5]], dtype=datatypes[5])
    attrs.coord2 = np.asarray(data[ll[5] :], dtype=datatypes[6])

    print(
        "attribute shapes: ",
        attrs.sdgbxindex.shape,
        attrs.xi.shape,
        attrs.radius.shape,
        attrs.msol.shape,
        attrs.coord3.shape,
        attrs.coord1.shape,
        attrs.coord2.shape,
    )

    return attrs


def totmass(radius, msol, RHO_L, RHO_SOL):
    """droplet totmass = mass of water + solute"""
    massconst = 4.0 / 3.0 * np.pi * radius * radius * radius * RHO_L
    density_factor = 1.0 - RHO_L / RHO_SOL
    totmass = msol * density_factor + massconst

    return totmass


def print_initsupers_infos(
    initsupers_filename, config_filename, constants_filename, grid_filename
):
    gbxvols = np.asarray(
        get_gbxvols_from_gridfile(
            grid_filename, constants_filename=constants_filename, isprint=False
        )
    )

    attrs = get_superdroplet_attributes(
        config_filename, constants_filename, initsupers_filename
    )

    xi = attrs.xi.flatten()
    vol = np.sum(gbxvols)
    numconc = np.sum(xi) / vol / 1e6  # [/cm^3]
    massconc = np.sum(attrs.msol.flatten() * xi) / vol * 1000  # [g m^-3]
    dropvol = 4 / 3 * np.pi * np.sum((attrs.radius.flatten() ** 3) * xi)
    m_w_conc = (
        dropvol * 1000 / vol * 1000
    )  # mass as if drops had density of water=1000Kg/m^3 [g m^3]

    inforstr = (
        "\n------ DOMAIN SUPERDROPLETS INFO ------\n"
        + "total droplet number conc: {:3g}".format(numconc)
        + " /cm^3\n"
        + "total droplet mass:        {:3g}".format(massconc)
        + " g/m^3\n"
        + "       as if water:        {:3g}".format(m_w_conc)
        + " g/m^3"
        + "\n------------------------------------\n"
    )
    print(inforstr)


def plot_initdistribs(attrs, gbxvols, gbxidxs):
    plt.rcParams.update({"font.size": 14})
    fig, axs = figure_setup(attrs.coord3, attrs.coord1, attrs.coord2)

    # create nbins evenly spaced in log10(r)
    nbins = 20
    minr, maxr = np.min(attrs.radius) / 10, np.max(attrs.radius) * 10
    hedgs = np.linspace(np.log10(minr), np.log10(maxr), nbins + 1)  # edges to lnr bins

    for idx in gbxidxs:
        vol = gbxvols[idx]
        sl = np.s_[attrs.sdgbxindex == idx]
        l0 = plot_radiusdistrib(axs[0], hedgs, attrs.radius[sl], attrs.xi[sl])
        l1 = plot_numconcdistrib(axs[1], hedgs, attrs.xi[sl], attrs.radius[sl], vol)
        l2 = plot_masssolutedistrib(
            axs[2], hedgs, attrs.xi[sl], attrs.radius[sl], attrs.msol[sl], vol
        )
        ls = plot_coorddistribs(axs, sl, hedgs, attrs)

    fig.tight_layout()

    return fig, axs, [l0, l1, l2, ls]


def plot_initGBxs_attrdistribs(
    config_filename,
    constants_filename,
    initsupers_filename,
    grid_filename,
    gbxs2plt,
    savefigpath,
    savefig,
    savelabel,
):
    """plot initial superdroplet distribution from initsupersfile binary
    of every gridbox with index in gbx2plts"""

    gbxvols = get_gbxvols_from_gridfile(
        grid_filename, constants_filename=constants_filename, isprint=False
    )
    attrs = get_superdroplet_attributes(
        config_filename, constants_filename, initsupers_filename
    )

    if isinstance(gbxs2plt, int):
        gbxidxs = [gbxs2plt]
        savename = "initGBx" + str(gbxs2plt) + "_distrib" + savelabel + ".png"
    elif gbxs2plt == "all":
        gbxidxs = np.unique(attrs.sdgbxindex)
        savename = "initallGBxs_distribs" + savelabel + ".png"
    else:
        gbxidxs = [int(g) for g in gbxs2plt]
        savename = "initGBxs_distribs" + savelabel + ".png"

    fig, axs, lines = plot_initdistribs(attrs, gbxvols, gbxidxs)

    fig.tight_layout()

    if savefig:
        savename = savefigpath / savename
        fig.savefig(savename, dpi=400, bbox_inches="tight", facecolor="w", format="png")
        print("Figure .png saved as: " + str(savename))


def figure_setup(coord3, coord1, coord2):
    ncoords = 3 - sum(not x.any() for x in [coord3, coord2, coord1])

    if ncoords == 0:
        fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(14, 4))
    elif ncoords == 1:
        fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(14, 8))
    else:
        fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(14, 8))

    fig.suptitle(str(ncoords) + "-D SDM Initial Superdroplet Conditions")
    axs = axs.flatten()
    if ncoords == 2:
        axs[-1].remove()

    return fig, axs


def log10r_frequency_distribution(radius, hedgs, wghts):
    """get distribution of data with weights 'wghts' against
    log10(r). Uses np.histogram to get frequency of a particular
    value of data that falls in each bin (with each bin defined
    by it's edges 'hedgs'). Return distirbution alongside the radius
    bin centers and widths in [m]"""

    if type(wghts) != np.ndarray:
        wghts = np.full(np.shape(radius), wghts)

    hist, hedgs = np.histogram(
        np.log10(radius), bins=hedgs, weights=wghts, density=None
    )

    # convert [m] to [micron]
    hedgs = (10 ** (hedgs)) * 1e6
    # radius bin widths [micron]
    hwdths = hedgs[1:] - hedgs[:-1]
    # radius bin centres [micron]
    hcens = (hedgs[1:] + hedgs[:-1]) / 2

    return hist, hedgs, hwdths, hcens


def plot_radiusdistrib(ax, hedgs, radius, xi):
    """get and plotthe superdroplet radius in each log10(r)
    bin and as a scatter on a twinx axis with their multiplicities"""

    l1 = ax.scatter(radius * 1e6, xi, zorder=1, label="multiplicities")

    ax2 = ax.twinx()
    hist, hedgs, hwdths, hcens = log10r_frequency_distribution(radius, hedgs, 1)
    l2 = ax2.step(
        hcens,
        hist,
        where="mid",
        alpha=0.8,
        zorder=0,
        color="grey",
        label="number distribution",
    )

    ax.set_xscale("log")
    ax.set_xlabel("radius, r, /\u03BCm")
    ax.set_yscale("log")

    ax.set_ylabel("superdroplet multiplicity")
    ax2.set_ylabel("superdroplet number distribution")

    if not ax.get_legend():
        ax.legend(loc="lower left")
        ax2.legend(loc="lower right")

    return [l1, l2]


def plot_numconcdistrib(ax, hedgs, xi, radius, vol):
    """get and plot frequency of real droplets in each log10(r) bin"""

    wghts = xi / vol / 1e6  # [cm^-3]
    hist, hedgs, hwdths, hcens = log10r_frequency_distribution(radius, hedgs, wghts)

    line = ax.step(hcens, hist, label="binned distribution", where="mid")
    ax.set_xscale("log")
    ax.set_xlabel("radius, r, /\u03BCm")
    ax.set_ylabel("real droplet number\nconcentration / cm$^{-3}$")

    if not ax.get_legend():
        ax.legend(loc="lower left")

    return line


def plot_masssolutedistrib(ax, hedgs, xi, radius, msol, vol):
    """get and plot frequency of real droplets in each log10(r) bin"""

    wghts = msol * xi / vol * 1000 / 1e6  # [g cm^-3]
    hist, hedgs, hwdths, hcens = log10r_frequency_distribution(radius, hedgs, wghts)

    line = ax.step(hcens, hist, where="mid")
    ax.set_xscale("log")
    ax.set_xlabel("radius, r, /\u03BCm")
    ax.set_ylabel("solute mass per unit volume / g cm$^{-3}$")

    return line


def plot_totmassdistrib(ax, hedgs, xi, radius, msol, vol, RHO_L, RHO_SOL):
    """get and plot frequency of real droplets in each log10(r) bin"""

    mass = totmass(radius, msol, RHO_L, RHO_SOL)
    wghts = mass * xi / vol * 1000 / 1e6  # [g cm^-3]
    hist, hedgs, hwdths, hcens = log10r_frequency_distribution(radius, hedgs, wghts)

    line = ax.step(hcens, hist, where="mid")
    ax.set_xscale("log")
    ax.set_xlabel("radius, r, /\u03BCm")
    ax.set_ylabel("droplet mass per unit volume / g cm$^{-3}$")

    return line


def scatter_totmass_solutemass(ax, radius, msol, RHO_L, RHO_SOL):
    mass = totmass(radius, msol, RHO_L, RHO_SOL)

    line = ax.scatter(mass * 1000, msol * 1000, marker="x")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("total mass / g")
    ax.set_ylabel("solute mass / g")

    return line


def plot_coorddistribs(axs, i2plt, hedgs, attrs):
    ls = []
    if attrs.coord3.any():
        ls.append(
            plot_coorddist(axs[3], hedgs, attrs.coord3[i2plt], attrs.radius[i2plt], 3)
        )
        if attrs.coord1.any():
            ls.append(
                plot_coorddist(
                    axs[4], hedgs, attrs.coord1[i2plt], attrs.radius[i2plt], 1
                )
            )
            if attrs.coord2.any():
                ls.append(
                    plot_coorddist(
                        axs[5], hedgs, attrs.coord2[i2plt], attrs.radius[i2plt], 2
                    )
                )
    return ls


def plot_coorddist(ax, hedgs, coord3, radius, coordnum):
    line = None
    if any(coord3):
        line = ax.scatter(radius * 1e6, coord3)

    ax.set_xscale("log")
    ax.set_xlabel("radius, r, /\u03BCm")
    ax.set_ylabel("superdroplet coord" + str(coordnum) + " / m")

    return line


def plot_initGBxs_dropletmasses(
    config_filename,
    constants_filename,
    initsupers_filename,
    grid_filename,
    gbxs2plt,
    savefigpath,
    savefig,
    savelabel,
):
    """plot initial superdroplet mass distributions
    from initsupersfile binary of every gridbox with index
    in gbx2plts"""

    gbxvols = get_gbxvols_from_gridfile(
        grid_filename, constants_filename=constants_filename, isprint=False
    )
    attrs = get_superdroplet_attributes(
        config_filename, constants_filename, initsupers_filename
    )
    inputs = initsupers_inputsdict(config_filename, constants_filename)

    if isinstance(gbxs2plt, int):
        gbxidxs = [gbxs2plt]
        savename = "initGBx" + str(gbxs2plt) + "_dropletmasses" + savelabel + ".png"
    elif gbxs2plt == "all":
        gbxidxs = np.unique(attrs.sdgbxindex)
        savename = "initallGBxs_dropletmasses" + savelabel + ".png"
    else:
        gbxidxs = [int(g) for g in gbxs2plt]
        savename = "initGBxs_dropletmasses" + savelabel + ".png"

    fig, axs, lines = plot_massdistribs(
        attrs, gbxvols, gbxidxs, inputs["RHO_L"], inputs["RHO_SOL"]
    )

    fig.tight_layout()

    if savefig:
        savename = savefigpath / savename
        fig.savefig(savename, dpi=400, bbox_inches="tight", facecolor="w", format="png")
        print("Figure .png saved as: " + str(savename))


def plot_massdistribs(attrs, gbxvols, gbxidxs, RHO_L, RHO_SOL):
    plt.rcParams.update({"font.size": 14})
    fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(14, 4))
    fig.suptitle("SDM Initial Superdroplet Mass Distributions")

    # create nbins evenly spaced in log10(r)
    nbins = 20
    minr, maxr = np.min(attrs.radius) / 10, np.max(attrs.radius) * 10
    hedgs = np.linspace(np.log10(minr), np.log10(maxr), nbins + 1)  # edges to lnr bins

    for idx in gbxidxs:
        vol = gbxvols[idx]
        sl = np.s_[attrs.sdgbxindex == idx]
        l0 = plot_totmassdistrib(
            axs[0],
            hedgs,
            attrs.xi[sl],
            attrs.radius[sl],
            attrs.msol[sl],
            vol,
            RHO_L,
            RHO_SOL,
        )
        l1 = plot_masssolutedistrib(
            axs[1], hedgs, attrs.xi[sl], attrs.radius[sl], attrs.msol[sl], vol
        )
        l2 = scatter_totmass_solutemass(
            axs[2], attrs.radius[sl], attrs.msol[sl], RHO_L, RHO_SOL
        )

    fig.tight_layout()

    return fig, axs, [l0, l1, l2]
