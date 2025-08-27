"""
Copyright (c) 2024 MPI-M, Clara Bayley


----- CLEO -----
File: create_initsuperdrops.py
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
from os.path import isfile
from .. import cxx2py, readconfigfile, writebinary
from ..gbxboundariesbinary_src.read_gbxboundaries import (
    read_dimless_gbxboundaries_binary,
)


class ManyAttrs:
    """store for lists of each attribute for all superdroplets"""

    def __init__(self):
        self.sdgbxindex = []
        self.xi = []
        self.radius = []
        self.msol = []
        self.coord3 = []
        self.coord1 = []
        self.coord2 = []

    def set_attrlists(self, a, b, c, d, e, f, g):
        self.sdgbxindex = a
        self.xi = b
        self.radius = c
        self.msol = d
        self.coord3 = e
        self.coord1 = f
        self.coord2 = g

    def extend_attrlists(self, mia):
        """use an instance of ManyAttrs (mia) to
        extend lists in this instance"""
        self.sdgbxindex.extend(mia.sdgbxindex)
        self.xi.extend(mia.xi)
        self.radius.extend(mia.radius)
        self.msol.extend(mia.msol)
        self.coord3.extend(mia.coord3)
        self.coord1.extend(mia.coord1)
        self.coord2.extend(mia.coord2)


def initsupers_inputsdict(config_filename, constants_filename):
    """create values from constants file & config file
    required as inputs to create initial
    superdroplet conditions"""

    consts = cxx2py.read_cxxconsts_into_floats(constants_filename)
    mconsts = cxx2py.derive_more_floats(consts)
    config = readconfigfile.read_configparams_into_floats(config_filename)

    inputs = {
        # for creating SD attribute distirbutions
        "nspacedims": int(config["nspacedims"]),
        "RHO_L": consts["RHO_L"],  # liquid density [Kg/m^3]
        "RHO_SOL": consts["RHO_SOL"],  # solute density [Kg/m^3]
        # for de-dimensionalising attributes
        "R0": consts["R0"],  # droplet radius lengthscale [m]
        "RHO0": mconsts["RHO0"],  # characteristic density scale [Kg/m^3]
        "MASS0": mconsts["MASS0"],  # characteristic mass scale [Kg]
        "COORD0": mconsts["COORD0"],  # z coordinate lengthscale [m]
    }

    return inputs


def is_sdgbxindex_correct(gridboxbounds, coord3, coord1, coord2, gbxindex, sdgbxindex):
    """rasises error if superdroplet coords [m] lie outside gridbox bounds [m]
    or if gridbox index not equal to superdroplet's."""

    errmsg = None
    for i, coord in enumerate([coord3, coord1, coord2]):
        if (coord < gridboxbounds[2 * i]).any():
            errmsg = (
                "superdroplet coord"
                + str(i + 1)
                + " below lower"
                + " limit of gridbox it's associated with"
            )
        elif (coord >= gridboxbounds[2 * i + 1]).any():
            errmsg = (
                "superdroplet coord"
                + str(i + 1)
                + " above upper"
                + " limit of gridbox it's associated with"
            )
    if errmsg:
        raise ValueError(errmsg)

    elif (sdgbxindex != gbxindex).any():
        errmsg = (
            "superdroplet gridbox index not the same as"
            + " gridbox it should be associated with"
        )
        raise ValueError(errmsg)


def dimless_superdropsattrs(
    nsupers,
    initattrsgen,
    inputs,
    gbxindex,
    gridboxbounds,
    NUMCONC,
    numconc_tolerance=0.0,
    isprint=False,
):
    """use superdroplet attribute generator "initattrsgen"
    and settings from config and consts files to
    make dimensionless superdroplet attributes"""

    # generate attributes
    sdgbxindex = [gbxindex] * nsupers
    xi, radius, msol = initattrsgen.generate_attributes(
        nsupers,
        inputs["RHO_SOL"],
        NUMCONC,
        gridboxbounds,
        numconc_tolerance=numconc_tolerance,
        isprint=isprint,
    )
    coord3, coord1, coord2 = initattrsgen.generate_coords(
        nsupers, inputs["nspacedims"], gridboxbounds
    )
    is_sdgbxindex_correct(gridboxbounds, coord3, coord1, coord2, gbxindex, sdgbxindex)

    # de-dimsionalise attributes
    radius = radius / inputs["R0"]
    msol = msol / inputs["MASS0"]
    coord3 = coord3 / inputs["COORD0"]
    coord1 = coord1 / inputs["COORD0"]
    coord2 = coord2 / inputs["COORD0"]

    attrs4gbx = ManyAttrs()
    attrs4gbx.set_attrlists(sdgbxindex, xi, radius, msol, coord3, coord1, coord2)

    return attrs4gbx


def create_allsuperdropattrs(
    nsupersdict,
    initattrsgen,
    gbxbounds,
    inputs,
    NUMCONC,
    numconc_tolerance=0.0,
    isprint=False,
):
    """returns lists for attributes of all SDs in domain called attrs"""

    attrs = ManyAttrs()  # lists of attrs for SDs in domain

    for gbxindex, gridboxbounds in gbxbounds.items():
        nsupers = nsupersdict[gbxindex]
        attrs4gbx = dimless_superdropsattrs(
            nsupers,
            initattrsgen,
            inputs,
            gbxindex,
            gridboxbounds,
            NUMCONC,
            numconc_tolerance=numconc_tolerance,
            isprint=isprint,
        )  # lists of attrs for SDs in gridbox

        attrs.extend_attrlists(attrs4gbx)

    return attrs


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

    return arr


def ctype_compatible_attrs(attrs):
    """make list from arrays of SD attributes that are
    compatible with c type expected by SDM e.g. unsigned
    long ints for xi, doubles for radius and msol"""

    datatypes = [np.uintc, np.uint, np.double, np.double]
    datatypes += [np.double] * 3  # coords datatype

    attrs.sdgbxindex = list(set_arraydtype(attrs.sdgbxindex, datatypes[0]))
    attrs.xi = list(set_arraydtype(attrs.xi, datatypes[1]))
    attrs.radius = list(set_arraydtype(attrs.radius, datatypes[2]))
    attrs.msol = list(set_arraydtype(attrs.msol, datatypes[3]))

    datalist = attrs.sdgbxindex + attrs.xi + attrs.radius + attrs.msol

    if any(attrs.coord3):
        # make coord3 compatible if there is data for it (>= 1-D model)
        attrs.coord3 = list(set_arraydtype(attrs.coord3, datatypes[4]))
        datalist += attrs.coord3

        if any(attrs.coord1):
            # make coord1 compatible if there is data for it and coord3 (>= 2-D model)
            attrs.coord1 = list(set_arraydtype(attrs.coord1, datatypes[4]))
            datalist += attrs.coord1

            if any(attrs.coord2):
                # make coord2 compatible if there is data for it and coord3 and coord2 (>= 3-D model)
                attrs.coord2 = list(set_arraydtype(attrs.coord2, datatypes[4]))
                datalist += attrs.coord2

    return datalist, datatypes


def check_datashape(data, ndata, nspacedims):
    """make sure each superdroplet attribute in data has length stated
    in ndata and that this length is compatible with the nummber of
    attributes and superdroplets expected given ndata"""

    err = ""
    if any([n != ndata[0] for n in ndata[: 4 + nspacedims]]):
        err += (
            "\n------ ERROR! -----\n"
            + "not all variables in data are same length, ndata = "
            + str(ndata[: 4 + nspacedims])
            + "\n---------------------\n"
        )

    if len(data) != np.sum(ndata):
        err += (
            "inconsistent dimensions of data: "
            + str(np.shape(data))
            + ", and"
            + " data per attribute: "
            + str(ndata)
            + ". data should be 1D with"
            + " shape: num_attributes * nsupers. data should be list of"
            + " [nsupers]*num_attributes."
        )

    if err:
        raise ValueError(err)


def nsupers_pergridboxdict(nsupers, gbxbounds):
    if isinstance(nsupers, int):
        nsupersdict = {}
        for key in gbxbounds.keys():
            nsupersdict[key] = nsupers
        return nsupersdict

    elif isinstance(nsupers, dict):
        if nsupers.keys() != gbxbounds.keys():
            errmsg = "keys for nsupers dict don't match gridbox indexes"
            raise ValueError(errmsg)
        else:
            return nsupers

    else:
        errmsg = "nsupers must be either dict or int"
        raise ValueError(errmsg)


def write_initsuperdrops_binary(
    initsupers_filename,
    initattrsgen,
    config_filename,
    constants_filename,
    grid_filename,
    nsupers,
    NUMCONC,
    numconc_tolerance=0.0,
    isprint=False,
):
    """de-dimensionalise attributes in initattrsgen and then write to
    to a binary file, "initsupers_filename", with some metadata"""

    if not isfile(grid_filename):
        errmsg = (
            "grid_filename not found, but must be"
            + " created before initsupersfile can be"
        )
        raise ValueError(errmsg)

    inputs = initsupers_inputsdict(config_filename, constants_filename)
    gbxbounds = read_dimless_gbxboundaries_binary(
        grid_filename, COORD0=inputs["COORD0"], isprint=False
    )
    nsupersdict = nsupers_pergridboxdict(nsupers, gbxbounds)

    attrs = create_allsuperdropattrs(
        nsupersdict,
        initattrsgen,
        gbxbounds,
        inputs,
        NUMCONC,
        numconc_tolerance=numconc_tolerance,
        isprint=isprint,
    )

    ndata = [
        len(dt)
        for dt in [
            attrs.sdgbxindex,
            attrs.xi,
            attrs.radius,
            attrs.msol,
            attrs.coord3,
            attrs.coord1,
            attrs.coord2,
        ]
    ]

    data, datatypes = ctype_compatible_attrs(attrs)
    check_datashape(data, ndata, inputs["nspacedims"])

    units = [b" ", b" ", b"m", b"g"]
    units += [b"m"] * 3  # coords units

    scale_factors = [1.0, 1.0, inputs["R0"], inputs["MASS0"]]
    scale_factors += [inputs["COORD0"]] * 3  # coords scale factors
    scale_factors = np.asarray(scale_factors, dtype=np.double)

    metastr = "Variables in this file are Superdroplet attributes:"
    if initattrsgen.coord3gen:
        if initattrsgen.coord1gen:
            if initattrsgen.coord2gen:
                metastr += " [sdgbxindex, xi, radius, msol, coord3, coord1, coord2]"
            else:
                metastr += " [sdgbxindex, xi, radius, msol, coord3, coord1]"
        else:
            metastr += " [sdgbxindex, xi, radius, msol, coord3]"
    else:
        metastr += " [sdgbxindex, xi, radius, msol]"

    writebinary.writebinary(
        initsupers_filename, data, ndata, datatypes, units, scale_factors, metastr
    )
