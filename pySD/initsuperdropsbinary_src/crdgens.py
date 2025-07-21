"""
----- CLEO -----
File: crdgens.py
Project: initsuperdropsbinary_src
Created Date: Friday 13th October 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
Copyright (c) 2023 MPI-M, Clara Bayley
-----
File Description:
various ways of generating spatial
coordinates for superdroplet initial
conditions
"""

import numpy as np

from ..gbxboundariesbinary_src import read_gbxboundaries as rgrid


def nsupers_at_domain_base(grid_filename, constants_filename, nsupers, zlim):
    """create dict for sd initialisation where nsupers
    only occur in gridboxes with upper bound <= zlim"""

    COORD0 = rgrid.get_COORD0_from_constsfile(constants_filename)
    gbxbounds, ndims = rgrid.read_dimless_gbxboundaries_binary(
        grid_filename, COORD0=COORD0, return_ndims=True, isprint=False
    )
    nsupersdict = {}
    for ii in gbxbounds.keys():
        gbx_zupper = gbxbounds[ii][1]  # z upper bound of gridbox
        if gbx_zupper <= zlim:
            nsupersdict[ii] = nsupers
        else:
            nsupersdict[ii] = 0

    return nsupersdict


def nsupers_at_domain_top(grid_filename, constants_filename, nsupers, zlim):
    """create dict for sd initialisation where nsupers
    only occur in gridboxes with lower bound >= zlim"""

    COORD0 = rgrid.get_COORD0_from_constsfile(constants_filename)
    gbxbounds, ndims = rgrid.read_dimless_gbxboundaries_binary(
        grid_filename, COORD0=COORD0, return_ndims=True, isprint=False
    )
    nsupersdict = {}
    for ii in gbxbounds.keys():
        gbx_zlower = gbxbounds[ii][0]  # z lower bound of gridbox
        if gbx_zlower >= zlim:
            nsupersdict[ii] = nsupers
        else:
            nsupersdict[ii] = 0

    return nsupersdict


class MonoCoordGen:
    """method to generate superdroplets with
    coord all equal to coord0"""

    def __init__(self, coord0):
        self.coord0 = coord0

    def __call__(self, nsupers, coordrange):
        """Returns coord for nsupers all
        with the value of coord0"""

        if self.coord0 >= coordrange[0] and self.coord0 < coordrange[1]:
            attrs = np.full(nsupers, self.coord0)
        else:
            attrs = np.array([])

        return attrs


class SampleCoordGen:
    """method to generate 'nsupers'
    no. of superdroplets' coord [m]
    by sampling in range bewteen coordspan"""

    def __init__(self, random):
        self.random = random

    def __call__(self, nsupers, coordrange):
        """Returns coord3 for nsupers
        sampled from coord3span [m]"""

        if not self.random:
            coord = np.linspace(coordrange[0], coordrange[1], nsupers, endpoint=False)
        else:
            coord = np.random.uniform(
                low=coordrange[0], high=coordrange[1], size=nsupers
            )

        return coord  # units [m]
