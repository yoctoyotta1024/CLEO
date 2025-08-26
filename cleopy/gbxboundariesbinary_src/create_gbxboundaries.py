"""
Copyright (c) 2024 MPI-M, Clara Bayley


----- CLEO -----
File: create_gbxboundaries.py
Project: gbxboundariesbinary_src
Created Date: Monday 16th October 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
"""

import numpy as np

from .. import cxx2py, writebinary


def get_COORD0_from_constsfile(constants_filename, returnconsts=False):
    """create values from constants file required as inputs to create initial
    superdroplet conditions"""

    consts = cxx2py.read_cxxconsts_into_floats(constants_filename)
    COORD0 = consts["TIME0"] * consts["W0"]

    if returnconsts:
        return COORD0, consts
    else:
        return COORD0


def linearlyspaced_halfcoords(grid):
    """returns linearly spaced half coords ie. 1D
    gridbox boundaries given 'grid' list of
    [min coord, max coord, delta coord] for x, y or z"""

    min, max, delta = grid

    return np.arange(min, max + delta, delta, dtype=np.double)


def get_dimless_halfcoords(grid, coord, COORD0):
    """given a list or array, return dimensionless halfcoords.
    For list call "halfcoords_from_coordlims" before
    de-dimensionalising. For array simply return array / COORD0"""

    if type(grid) not in [list, np.ndarray]:
        raise ValueError("input " + coord + " grid is neither list or array ")

    elif isinstance(grid, list):
        grid = linearlyspaced_halfcoords(grid)

    elif isinstance(grid, np.ndarray):
        grid = np.array(grid, dtype=np.double)

    return grid / COORD0


def check_halfcoords(grid, coord):
    """check that x , y or z grid limits are for at least
    1 cell with strictly monotonically increasing bounds"""

    if coord not in ["z", "x", "y"]:
        raise ValueError("coord should be x, y or z")

    criteria = len(grid) >= 2 and np.all(np.diff(grid) > 0)
    if criteria is False:
        errmsg = str(grid) + " does not meet criteria for " + coord + " halfcoords"
        raise ValueError(errmsg)


def gridboxboundaries_from_halfcoords(allhalfcoords, ngridboxes):
    """returns gbxbounds dictionary. Each key of gbxbounds is a gridbox's
    index and the corresponding value is the gridbox's
    [zmin, zmax, xmin, xmax, ymin, ymax] boundaries"""

    zhalfs = allhalfcoords[0]
    xhalfs = allhalfcoords[1]
    yhalfs = allhalfcoords[2]

    gbxbounds, ii = {}, 0
    for j in range(len(yhalfs) - 1):
        for i in range(len(xhalfs) - 1):
            for k in range(len(zhalfs) - 1):
                zbounds = [zhalfs[k], zhalfs[k + 1]]
                xbounds = [xhalfs[i], xhalfs[i + 1]]
                ybounds = [yhalfs[j], yhalfs[j + 1]]
                gbxbounds[ii] = zbounds + xbounds + ybounds
                ii += 1

    if len(gbxbounds) != ngridboxes:
        raise ValueError("not enough gbxbounds for ngridboxes")
    else:
        print("created boundaries for", ngridboxes, "gridboxes")

    return gbxbounds


def dimless_gridboxboundaries(zgrid, xgrid, ygrid, COORD0):
    """use zgrid, xgrid and ygrid lists or arrays to create half coords
    of gridboxes (ie. their boundaries). Return single list of zhalf,
    then xhalf, then yhalf coords and a list of len(3) which states how
    many of each z, x and y coords there are"""

    allhalfcoords, ndims = [], []  # z, x and y half coords in 1 continuous list
    for grid, coord in zip([zgrid, xgrid, ygrid], ["z", "x", "y"]):
        halfcoords = get_dimless_halfcoords(grid, coord, COORD0)
        check_halfcoords(halfcoords, coord)

        allhalfcoords.append(halfcoords)
        ndims.append(len(halfcoords) - 1)

    gbxbounds = gridboxboundaries_from_halfcoords(allhalfcoords, np.prod(ndims))

    ndims = np.asarray(ndims, dtype=np.uint)
    gbxindicies = np.array(list(gbxbounds.keys()), dtype=np.uintc).flatten()
    gbxboundsdata = np.array(list(gbxbounds.values()), dtype=np.double).flatten()

    return ndims, gbxindicies, gbxboundsdata


def set_arraydtype(arr, dtype):
    og = type(arr[0])
    if og != dtype:
        arr = np.array(arr, dtype=dtype)

        warning = (
            "WARNING! dtype of attributes is being changed!"
            + " from "
            + str(og)
            + " to "
            + str(dtype)
        )
        raise ValueError(warning)

    return list(arr)


def ctype_compatible_gridboxboundaries(ndims, idxs, bounds):
    """check type of gridbox boundaries data is compatible
    with c type double. If not, change type and raise error"""

    datatypes = [np.uint, np.uintc, np.double]

    ndims = set_arraydtype(ndims, datatypes[0])
    idxs = set_arraydtype(idxs, datatypes[1])
    bounds = set_arraydtype(bounds, datatypes[2])

    datalist = ndims + idxs + bounds

    return datalist, datatypes


def write_gridboxboundaries_binary(
    grid_filename, zgrid, xgrid, ygrid, constants_filename
):
    """zgrid, xgrid and ygrid can either be list of
    [min coord, max coord, delta] or they can be arrays of
    their half coordinates (ie. gridbox boundaries). If the former,
    first create half coordinates. In both cases, de-dimensionalise
    and write the boundaries to a binary file, "filename" """

    COORD0 = get_COORD0_from_constsfile(constants_filename)

    ndims, gbxindicies, gbxboundsdata = dimless_gridboxboundaries(
        zgrid, xgrid, ygrid, COORD0
    )
    zxy = [len(zgrid) - 1, len(xgrid) - 1, len(ygrid) - 1]
    zxy = [str(x) for x in zxy]
    metastr = (
        "Variables in this file are ndims in (z,x,y), then the "
        + str(np.prod(ndims))
        + " gridbox indicies followed by the"
        + " [zmin, zmax, xmin, xmax, ymin, ymax]"
        + " coordinates for each gridbox's boundaries."
        + " Grid has dimensions "
        + zxy[0]
        + "x"
        + zxy[1]
        + "x"
        + zxy[2]
    )

    ndata = [len(dt) for dt in [ndims, gbxindicies, gbxboundsdata]]
    data, datatypes = ctype_compatible_gridboxboundaries(
        ndims, gbxindicies, gbxboundsdata
    )
    scale_factors = np.array([1.0, 1.0, COORD0], dtype=np.double)
    units = [b" ", b" ", b"m"]

    writebinary.writebinary(
        grid_filename, data, ndata, datatypes, units, scale_factors, metastr
    )
