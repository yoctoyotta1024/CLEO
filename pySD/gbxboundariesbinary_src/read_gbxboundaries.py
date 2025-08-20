"""
Copyright (c) 2024 MPI-M, Clara Bayley


----- CLEO -----
File: read_gbxboundaries.py
Project: gbxboundariesbinary_src
Created Date: Wednesday 17th January 2024
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
from pathlib import Path

from .create_gbxboundaries import get_COORD0_from_constsfile
from ..readbinary import readbinary


def get_gridboxboundaries(
    grid_filename, COORD0=False, constants_filename="", isprint=True
):
    """get gridbox boundaries from binary file and
    re-dimensionalise usign COORD0 const from constants_filename"""

    if not COORD0:
        COORD0 = get_COORD0_from_constsfile(constants_filename)

    gbxbounds = read_dimless_gbxboundaries_binary(
        grid_filename, COORD0, isprint=isprint
    )

    zhalf, xhalf, yhalf = halfcoords_from_gbxbounds(gbxbounds, isprint=isprint)

    return zhalf, xhalf, yhalf


def get_domainvol_from_grid_filename(
    grid_filename, COORD0=False, constants_filename=""
):
    """get total domain volume from binary file"""

    zhalf, xhalf, yhalf = get_gridboxboundaries(
        grid_filename, COORD0=COORD0, constants_filename=constants_filename
    )

    return calc_domainvol(zhalf, xhalf, yhalf)


def get_gbxvols_from_gridfile(
    grid_filename, COORD0=False, constants_filename="", isprint=True
):
    """get total domain volume from binary file"""

    if not COORD0:
        COORD0 = get_COORD0_from_constsfile(constants_filename)

    gbxbounds = read_dimless_gbxboundaries_binary(
        grid_filename, COORD0, isprint=isprint
    )

    return calc_gridboxvols(gbxbounds)


def fullcell(halfcell):
    """return fullcell corods (cell centres)
    given half coords in a direction"""

    return (halfcell[1:] + halfcell[:-1]) / 2


def cellwidth(halfcell):
    """return cell spacing (width) given half coords
    in a direction"""

    delta = abs(halfcell[1:] - halfcell[:-1])

    return delta


def fullcell_fromhalfcoords(zhalf, xhalf, yhalf):
    zfull = fullcell(zhalf)
    xfull = fullcell(xhalf)
    yfull = fullcell(yhalf)

    return zfull, xfull, yfull


def fullcoords_forallgridboxes(gbxbounds, ndims):
    """returns (x,y,z) centres of gridboxes in domain"""

    zhalf, xhalf, yhalf = halfcoords_from_gbxbounds(gbxbounds, isprint=False)
    zfull, xfull, yfull = fullcell_fromhalfcoords(zhalf, xhalf, yhalf)

    zfullcoords = np.tile(
        zfull, int(ndims[1] * ndims[2])
    )  # zfull of every gridbox in order of gbxindex
    xfullcoords = np.tile(np.repeat(xfull, ndims[0]), int(ndims[2]))
    yfullcoords = np.repeat(yfull, ndims[0] * ndims[1])

    return zfullcoords, xfullcoords, yfullcoords


def coords_forgridboxfaces(gbxbounds, ndims, face):
    """returns (x,y,z) coordinates of gridboxes faces
    in a particular direction"""

    ndims = [int(n) for n in ndims]
    zhalf, xhalf, yhalf = halfcoords_from_gbxbounds(gbxbounds, isprint=False)
    zfull, xfull, yfull = fullcell_fromhalfcoords(zhalf, xhalf, yhalf)

    if face == "z":
        nz, nx, ny = ndims[0] + 1, ndims[1], ndims[2]
        zfaces = np.tile(zhalf, nx * ny)  # zfull of every gridbox in order of gbxindex
        xfulls = np.tile(np.repeat(xfull, nz), ny)
        yfulls = np.repeat(yfull, nz * nx)
        return zfaces, xfulls, yfulls

    elif face == "x":
        nz, nx, ny = ndims[0], ndims[1] + 1, ndims[2]
        zfulls = np.tile(zfull, nx * ny)  # zfull of every gridbox in order of gbxindex
        xfaces = np.tile(np.repeat(xhalf, nz), ny)
        yfulls = np.repeat(yfull, nz * nx)
        return zfulls, xfaces, yfulls

    elif face == "y":
        nz, nx, ny = ndims[0], ndims[1], ndims[2] + 1
        zfulls = np.tile(zfull, nx * ny)  # zfull of every gridbox in order of gbxindex
        xfulls = np.tile(np.repeat(xfull, nz), ny)
        yfaces = np.repeat(yhalf, nz * nx)
        return zfulls, xfulls, yfaces


def read_dimless_gbxboundaries_binary(
    filename, COORD0=False, return_ndims=False, isprint=True
):
    """return dictionary for gbx indicies to gbx boundaries by
    reading binary file. Return dimensionless version if COORD0
    not give (=False)."""

    data, ndata_pervar = readbinary(filename, isprint=isprint)
    datatypes = [np.uint, np.uintc, np.double]

    ll = [0, 0]  # indexs for division of data list between each variable
    for n in range(1, len(ndata_pervar)):
        ll[n - 1] = np.sum(ndata_pervar[:n])

    ndims = np.asarray(data[: ll[0]], dtype=datatypes[0])
    gbxidxs = np.asarray(data[ll[0] : ll[1]], dtype=datatypes[1])

    ngridboxes = int(np.prod(ndims))
    if len(gbxidxs) != ngridboxes:
        err = "number of gridbox indexes not consistent with (z,x,y) dims"
        raise ValueError(err)

    boundsdata = np.asarray(data[ll[1] :], dtype=datatypes[2])
    boundsdata = np.reshape(boundsdata, [ngridboxes, len(boundsdata) // ngridboxes])

    if COORD0:
        boundsdata = boundsdata * COORD0

    gbxbounds = {gbxidxs[i]: boundsdata[i] for i in range(ngridboxes)}

    if return_ndims:
        return gbxbounds, ndims
    else:
        return gbxbounds


def halfcoords_from_gbxbounds(gbxbounds, isprint=True):
    """returns half coords of gbx boundaries in lists obtained
    from gbxbounds dictionary"""

    boundsdata = np.asarray(list(gbxbounds.values()))

    zhalf = np.unique(np.sort(boundsdata[:, 0]))
    zhalf = np.append(zhalf, np.amax(boundsdata[:, 1]))

    xhalf = np.unique(np.sort(boundsdata[:, 2]))
    xhalf = np.append(xhalf, np.amax(boundsdata[:, 3]))

    yhalf = np.unique(np.sort(boundsdata[:, 4]))
    yhalf = np.append(yhalf, np.amax(boundsdata[:, 5]))

    if isprint:
        print("zhalf: ", zhalf)
        print("xhalf: ", xhalf)
        print("yhalf: ", yhalf)

    return zhalf, xhalf, yhalf


def plot_gridboxboundaries(
    constants_filename,
    grid_filename,
    savefig=False,
    savefigpath=None,
    savelabel="",
):
    plt.rcParams.update({"font.size": 14})

    zhalf, xhalf, yhalf = get_gridboxboundaries(
        grid_filename, constants_filename=constants_filename
    )

    halfs = [zhalf, xhalf, yhalf]
    fulls, deltas = [], []
    for half in halfs:
        fulls.append(fullcell(half))
        deltas.append(cellwidth(half))

    fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(10, 5))

    for i, crd in enumerate(["z", "x", "y"]):
        xlims = [0.8 * np.amin(deltas[i]) / 1000, np.amax(deltas[i]) / 1000 * 1.2]
        ylims = [np.amin(halfs[i]) / 1000, np.amax(halfs[i]) / 1000]
        axs[i].scatter(deltas[i] / 1000, fulls[i] / 1000, color="k", label="centres")
        axs[i].hlines(
            halfs[i] / 1000,
            xlims[0],
            xlims[1],
            color="grey",
            alpha=0.8,
            linewidth=0.8,
            label="boundaries",
        )
        axs[i].set_xlim(xlims)
        axs[i].set_ylim(ylims)
        axs[i].set_xlabel("gridbox spacing, \u0394 " + crd + " /km")
        axs[i].set_ylabel("gridbox centres, " + crd + "f /km")
        axs[i].legend()

    fig.tight_layout()

    if savefig:
        savename = savefigpath / Path(f"gridboxboundaries{savelabel}.png")
        fig.savefig(
            savename,
            dpi=400,
            bbox_inches="tight",
            facecolor="w",
            format="png",
        )
        print("Figure .png saved as: " + str(savename))


def calc_domainvol(zhalf, xhalf, yhalf):
    widths = []
    for half in [zhalf, xhalf, yhalf]:
        widths.append(np.amax(half) - np.amin(half))

    domainvol = np.prod(widths)

    return domainvol


def calc_gridboxvols(gbxbounds):
    gbxvols = []
    for gbxindex, bounds in gbxbounds.items():
        zwidth = bounds[1] - bounds[0]
        xwidth = bounds[3] - bounds[2]
        ywidth = bounds[5] - bounds[4]

        gbxvols.append(zwidth * xwidth * ywidth)

    return gbxvols


def domaininfo(gbxbounds, isprint=True):
    zhalf, xhalf, yhalf = halfcoords_from_gbxbounds(gbxbounds, isprint=isprint)
    domainvol = calc_domainvol(zhalf, xhalf, yhalf)

    gridboxvols = calc_gridboxvols(gbxbounds)
    ngridboxes = len(gridboxvols)

    return domainvol, gridboxvols, ngridboxes


def grid_dimensions(gbxbounds):
    zhalf, xhalf, yhalf = halfcoords_from_gbxbounds(gbxbounds, isprint=False)

    if len(zhalf) == 2:
        griddims = 0
    elif len(xhalf) == 2:
        griddims = 1
    elif len(yhalf) == 2:
        griddims = 2
    else:
        griddims = 3

    extents, spacings = [], []
    for half in [zhalf, xhalf, yhalf]:
        extents.append([np.amin(half), np.amax(half)])
        spacings.append(abs(half[1:] - half[:-1]))

    return extents, spacings, griddims


def print_domain_info(constants_filename, grid_filename):
    """prints information about domain read from grid_filename and constants file"""

    isprint = True
    COORD0 = get_COORD0_from_constsfile(constants_filename)
    gbxbounds = read_dimless_gbxboundaries_binary(
        grid_filename, COORD0, isprint=isprint
    )

    domainvol, gridboxvols, ngridboxes = domaininfo(gbxbounds, isprint=isprint)
    xtns, spacings, griddims = grid_dimensions(gbxbounds)
    ztot = abs(xtns[0][0] - xtns[0][1])
    xtot = abs(xtns[1][0] - xtns[1][1])
    ytot = abs(xtns[2][0] - xtns[2][1])

    inforstr = (
        "\n------ DOMAIN / GRIDBOXES INFO ------\n"
        + "------------- "
        + str(griddims)
        + "-D MODEL -------------\n"
        + "domain dimensions: ({:3g}x{:3g}x{:3g})m^3\n".format(ztot, xtot, ytot)
        + "domain no. gridboxes: "
        + str(len(spacings[0]))
        + "x"
        + str(len(spacings[1]))
        + "x"
        + str(len(spacings[2]))
        + "\n"
        + "domain z limits: ({:3g},{:3g})m\n".format(np.amin(xtns[0]), np.amax(xtns[0]))
        + "domain x limits: ({:3g}, {:3g})m\n".format(
            np.amin(xtns[1]), np.amax(xtns[1])
        )
        + "domain y limits: ({:3g}, {:3g})m\n".format(
            np.amin(xtns[2]), np.amax(xtns[2])
        )
        + "mean gridbox z spacing: {:3g} m\n".format(np.mean(spacings[0]))
        + "mean gridbox x spacing: {:3g} m\n".format(np.mean(spacings[1]))
        + "mean gridbox y spacing: {:3g} m\n".format(np.mean(spacings[2]))
        + "mean gridbox volume: {:3g}".format(np.mean(gridboxvols))
        + " m^3\n"
        + "total domain volume: {:3g} m^3\n".format(domainvol)
        + "total no. gridboxes: "
        + str(ngridboxes)
        + "\n------------------------------------\n"
    )
    print(inforstr)
