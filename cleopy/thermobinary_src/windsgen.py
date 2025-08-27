"""
Copyright (c) 2024 MPI-M, Clara Bayley


----- CLEO -----
File: windsgen.py
Project: thermobinary_src
Created Date: Wednesday 26th March 2025
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
Various ways of generating wind fields (wvel, vvel and uvel) to read into CLEO
"""

import numpy as np
from .create_thermodynamics import thermoinputsdict
from ..gbxboundariesbinary_src import read_gbxboundaries as rgrid


def empty_winds():
    WINDSDATA = {
        "WVEL": np.array([]),
        "UVEL": np.array([]),
        "VVEL": np.array([]),
    }

    return WINDSDATA


def constant_winds(ndims, ntime, WVEL, UVEL, VVEL):
    """create dictionary for winds, array are
    empty by default, or given constant values WVEl, UVEL and VVEL.
    Here, shape_[X]face = no. data for wind velocity component
    defined on gridbox [X] faces"""

    WINDSDATA = empty_winds()

    if WVEL is not None:
        shape_zface = int((ndims[0] + 1) * ndims[1] * ndims[2] * ntime)
        WINDSDATA["WVEL"] = np.full(shape_zface, WVEL)

        if UVEL is not None:
            shape_xface = int((ndims[1] + 1) * ndims[2] * ndims[0] * ntime)
            WINDSDATA["UVEL"] = np.full(shape_xface, UVEL)

            if VVEL is not None:
                shape_yface = int((ndims[2] + 1) * ndims[0] * ndims[1] * ntime)
                WINDSDATA["VVEL"] = np.full(shape_yface, VVEL)

    return WINDSDATA


def divfree_flowfield2D(
    wmax, zlength, xlength, rhotilda_zfaces, rhotilda_xfaces, gbxbounds, ndims
):
    zfaces, xcens_z = rgrid.coords_forgridboxfaces(gbxbounds, ndims, "z")[0:2]
    zcens_x, xfaces = rgrid.coords_forgridboxfaces(gbxbounds, ndims, "x")[0:2]

    ztilda = zlength / np.pi
    xtilda = xlength / (2 * np.pi)
    wamp = 2 * wmax

    WVEL = wamp / rhotilda_zfaces
    WVEL = WVEL * np.sin(zfaces / ztilda) * np.sin(xcens_z / xtilda)

    UVEL = wamp / rhotilda_xfaces * xtilda / ztilda
    UVEL = UVEL * np.cos(zcens_x / ztilda) * np.cos(xfaces / xtilda)

    return WVEL, UVEL


class ConstUniformWinds:
    """create winds that's constant in time and uniform throughout the domain"""

    def __init__(
        self,
        WVEL,
        UVEL,
        VVEL,
    ):
        self.WVEL = WVEL  # vertical velocity [m/s]
        self.UVEL = UVEL  # horizontal eastwards velocity [m/s]
        self.VVEL = VVEL  # horizontal northwards velocity [m/s]

    def generate_winds(self, gbxbounds, ndims, ntime, THERMODATA):
        return constant_winds(ndims, ntime, self.WVEL, self.UVEL, self.VVEL)


class Simple2DFlowField:
    """create winds that are constant in time in a 2D (z,x) dependent flow field
    assumign dry hydrostatic adiabat for density profile.
    Equations derived from Arabas et al. 2015 (sect 2.1)"""

    def __init__(
        self,
        config_filename,
        constants_filename,
        WMAX,
        Zlength,
        Xlength,
        VVEL,
    ):
        inputs = thermoinputsdict(config_filename, constants_filename)

        self.WMAX = WMAX  # max velocities constant
        self.Zlength = Zlength  # wavelength of velocity modulation in z direction [m]
        self.Xlength = Xlength  # wavelength of velocity modulation in x direction [m]
        self.VVEL = VVEL  # horizontal (y) velocity

        self.RGAS_DRY = inputs["RGAS_DRY"]
        self.RGAS_V = inputs["RGAS_V"]
        self.RHO0 = inputs["RHO0"]

    def wvel_uvel_from_flowfield(self, THERMODATA, gbxbounds, ndims):
        PRESS, TEMP = THERMODATA["PRESS"][0], THERMODATA["TEMP"][0]
        qvap = THERMODATA["qvap"][0]
        rho_dry = PRESS / (TEMP * (self.RGAS_DRY + qvap * self.RGAS_V))
        rhotilda = rho_dry / self.RHO0  # scalar value is the same over entire domain

        WVEL, UVEL = divfree_flowfield2D(
            self.WMAX, self.Zlength, self.Xlength, rhotilda, rhotilda, gbxbounds, ndims
        )
        return WVEL, UVEL

    def generate_winds(self, gbxbounds, ndims, ntime, THERMODATA):
        WINDSDATA = empty_winds()

        if self.WMAX is not None:
            WVEL, UVEL = self.wvel_uvel_from_flowfield(THERMODATA, gbxbounds, ndims)
            WINDSDATA["WVEL"] = np.tile(WVEL, ntime)
            WINDSDATA["UVEL"] = np.tile(UVEL, ntime)

            if self.VVEL is not None:
                shape_yface = int((ndims[2] + 1) * ndims[0] * ndims[1] * ntime)
                WINDSDATA["VVEL"] = np.full(shape_yface, self.VVEL)

        return WINDSDATA


class DryAdiabat2DFlowField:
    """create winds that are constant in time in a 2D (z,x) dependent flow field.
    Divergence free field with density dependence assuming dry adiabat.
    Equations derived from Arabas et al. 2015 (sect 2.1)"""

    def __init__(
        self,
        WMAX,
        Zlength,
        Xlength,
        VVEL,
        dryadiabat,
    ):
        self.WMAX = WMAX  # max velocities constant
        self.Zlength = Zlength  # wavelength of velocity modulation in z direction [m]
        self.Xlength = Xlength  # wavelength of velocity modulation in x direction [m]
        self.VVEL = VVEL  # horizontal (y) velocity

        self.dryadiabat = dryadiabat

    def wvel_uvel_from_flowfield(self, gbxbounds, ndims):
        zfaces = rgrid.coords_forgridboxfaces(gbxbounds, ndims, "z")[0]
        xfaces = rgrid.coords_forgridboxfaces(gbxbounds, ndims, "x")[1]
        rhotilda_zfaces = self.dryadiabat.rhotilda(zfaces)
        rhotilda_xfaces = self.dryadiabat.rhotilda(xfaces)

        WVEL, UVEL = divfree_flowfield2D(
            self.WMAX,
            self.Zlength,
            self.Xlength,
            rhotilda_zfaces,
            rhotilda_xfaces,
            gbxbounds,
            ndims,
        )
        return WVEL, UVEL

    def generate_winds(self, gbxbounds, ndims, ntime, THERMODATA):
        WINDSDATA = empty_winds()

        if self.WMAX is not None:
            WVEL, UVEL = self.wvel_uvel_from_flowfield(gbxbounds, ndims)
            WINDSDATA["WVEL"] = np.tile(WVEL, ntime)
            WINDSDATA["UVEL"] = np.tile(UVEL, ntime)

            if self.VVEL is not None:
                shape_yface = int((ndims[2] + 1) * ndims[0] * ndims[1] * ntime)
                WINDSDATA["VVEL"] = np.full(shape_yface, self.VVEL)

        return WINDSDATA


class SinusoidalUpdraught:
    """creates constant in time wind field with sinusoidal vertial velocity profile 'wwel'
    and optionally uniform horizontal (uvel, vvel) wind field"""

    def __init__(
        self,
        WMAX,
        UVEL,
        VVEL,
        Wlength,
    ):
        self.WMAX = WMAX  # vertical (coord3) velocity [m/s]
        self.UVEL = UVEL  # horizontal eastwards (coord1) velocity [m/s]
        self.VVEL = VVEL  # horizontal northwards (coord2) velocity [m/s]
        self.Wlength = Wlength  # [m] use constant W (Wlength=0.0), or sinusoidal 1-D profile below cloud base

    def wvel_profile(self, gbxbounds, ndims, ntime):
        """returns updraught (w always >=0.0) sinusoidal
        profile with amplitude WMAX and wavelength 2*Wlength"""

        zfaces = rgrid.coords_forgridboxfaces(gbxbounds, ndims, "z")[0]
        WVEL = self.WMAX * np.sin(np.pi * zfaces / (2 * self.Wlength))

        WVEL[WVEL < 0.0] = 0.0

        return np.tile(WVEL, ntime)

    def generate_winds(self, gbxbounds, ndims, ntime, THERMODATA):
        WINDSDATA = constant_winds(ndims, ntime, self.WMAX, self.UVEL, self.VVEL)
        if self.Wlength > 0.0:
            WINDSDATA["WVEL"] = self.wvel_profile(gbxbounds, ndims, ntime)

        return WINDSDATA
