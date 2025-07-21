"""
----- CLEO -----
File: pygbxsdat.py
Project: sdmout_src
Created Date: Tuesday 24th October 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
Copyright (c) 2023 MPI-M, Clara Bayley
-----
File Description:
for reading data from gridbox (gbx)
boudnaries binary inputs files
"""

import numpy as np
from ..gbxboundariesbinary_src import read_gbxboundaries as rgrid
from ..cxx2py import print_dict_statement


def get_gridboxes(grid_filename, COORD0, isprint=True):
    gbxbounds, ndims = rgrid.read_dimless_gbxboundaries_binary(
        grid_filename, COORD0=COORD0, return_ndims=True, isprint=False
    )
    zhalf, xhalf, yhalf = rgrid.halfcoords_from_gbxbounds(gbxbounds, isprint=False)
    domainvol, gbxvols, ngrid = rgrid.domaininfo(gbxbounds, isprint=False)

    gbxs = {
        "ngrid": ngrid,  # number of gridboxes
        "ndims": np.flip(ndims),  # dimensions (no. gridboxes in [y,x,z] direction)
        "domainvol": domainvol,
        "domainarea": domainvol
        / (np.amax(zhalf) - np.amin(zhalf)),  # x-y plane horizontal are
        "gbxvols": gbxvols,  # list of vols of each gbx
        "zhalf": zhalf,  # half cell coords (boundaries)
        "zfull": rgrid.fullcell(zhalf),  # full cell coords (centres)
        "xhalf": xhalf,  # half cell coords (boundaries)
        "xfull": rgrid.fullcell(xhalf),  # full cell coords (centres)
        "yhalf": yhalf,  # half cell coords (boundaries)
        "yfull": rgrid.fullcell(yhalf),  # full cell coords (centres)
    }
    gbxs["gbxvols"] = np.reshape(gbxs["gbxvols"], gbxs["ndims"])

    xxh, zzh = np.meshgrid(
        gbxs["xhalf"], gbxs["zhalf"], indexing="ij"
    )  # dims [xdims, zdims] [m]
    xxf, zzf = np.meshgrid(
        gbxs["xfull"], gbxs["zfull"], indexing="ij"
    )  # dims [xdims, zdims] [m]
    gbxs["xxh"], gbxs["zzh"] = xxh, zzh
    gbxs["xxf"], gbxs["zzf"] = xxf, zzf

    if isprint:
        print_dict_statement(grid_filename, "gbxs", gbxs)

    return gbxs
