"""
Copyright (c) 2024 MPI-M, Clara Bayley


----- CLEO -----
File: thermogen.py
Project: thermobinary_src
Created Date: Monday 16th October 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
Last Modified: Tuesday 7th May 2024
Modified By: CB
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
"""

import numpy as np
from scipy import integrate
from .. import cxx2py
from .create_thermodynamics import thermoinputsdict
from ..gbxboundariesbinary_src import read_gbxboundaries as rgrid


def get_Mrratio_from_constants_filename(constants_filename):
    consts = cxx2py.read_cxxconsts_into_floats(constants_filename)
    mconsts = cxx2py.derive_more_floats(consts)

    return mconsts["Mr_ratio"]


def saturation_press(TEMP):
    """Calculate the equilibrium vapor pressure of water over
    liquid water ie. the saturation pressure (psat) given the
    temperature [K]. Equation from Bjorn Steven's "make_tetens"
    function in module "moist_thermodynamics.saturation_vapour_pressures"
    available on gitlab. Original paper "Murray, F. W. On the
    Computation of Saturation Vapor Pressure. Journal of Applied
    Meteorology and Climatology 6, 203–204 (1967)."""

    Aconst = 17.4146
    Bconst = 33.639
    TREF = 273.16  # Triple point temperature [K] of water
    PREF = 611.655  # Triple point pressure [Pa] of water

    if np.any(TEMP <= 0.0):
        raise ValueError("psat ERROR: T must be larger than 0K." + " T = " + str(TEMP))

    return PREF * np.exp(Aconst * (TEMP - TREF) / (TEMP - Bconst))  # [Pa]


def relh2qvap(press, temp, relh, Mr_ratio):
    """convert relative humidity [%] (relh) into vapour mass
    mixing ratio (qvap) given ambient temperature and pressure
    and ratio of molecular masses: vapour/air"""

    vapourpress = saturation_press(temp) * relh / 100.0  # [Pa]

    qvap = Mr_ratio * vapourpress / (press - vapourpress)  # dimensionless [Kg/kg]

    return qvap


def sratio2qvap(sratio, press, temp, Mr_ratio):
    psat = saturation_press(temp)

    qvap = Mr_ratio * sratio
    qvap = qvap / (press / psat - 1)

    return qvap


def qparams_to_qvap(method, params, Mr_ratio, PRESS, TEMP):
    """returns qvaps given list of qvaps, supersaturation ratios
    or relative humidities"""

    if method == "qvap":
        qparams = params
        return qparams

    elif method == "sratio":
        qparams = []
        for sratio in params:
            qparams.append(sratio2qvap(sratio, PRESS, TEMP, Mr_ratio))
        return qparams

    elif method == "relh":
        qparams = []
        for relh in params:
            qparams.append(relh2qvap(PRESS, TEMP, relh, Mr_ratio))
        return qparams

    else:
        raise ValueError("valid method not given to generate qvap")


def constant_winds(ndims, ntime, THERMODATA, WVEL, UVEL, VVEL):
    """add arrays to thermodata dictionary for winds, array are
    empty by default, or given constant values WVEl, UVEL and VVEL.
    Here, shape_[X]face = no. data for wind velocity component
    defined on gridbox [X] faces"""

    for VEL in ["WVEL", "UVEL", "VVEL"]:
        THERMODATA[VEL] = np.array([])

    if WVEL is not None:
        shape_zface = int((ndims[0] + 1) * ndims[1] * ndims[2] * ntime)
        THERMODATA["WVEL"] = np.full(shape_zface, WVEL)

        if UVEL is not None:
            shape_xface = int((ndims[1] + 1) * ndims[2] * ndims[0] * ntime)
            THERMODATA["UVEL"] = np.full(shape_xface, UVEL)

            if VVEL is not None:
                shape_yface = int((ndims[2] + 1) * ndims[0] * ndims[1] * ntime)
                THERMODATA["VVEL"] = np.full(shape_yface, VVEL)

    return THERMODATA


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


class ConstUniformThermo:
    """create thermodynamics that's constant in
    time and uniform throughout the domain"""

    def __init__(
        self,
        PRESS,
        TEMP,
        qvap,
        qcond,
        WVEL,
        UVEL,
        VVEL,
        relh=False,
        constants_filename="",
    ):
        self.PRESS = PRESS  # pressure [Pa]
        self.TEMP = TEMP  # temperature [T]

        if relh:
            Mr_ratio = get_Mrratio_from_constants_filename(constants_filename)
            self.qvap = relh2qvap(
                PRESS, TEMP, relh, Mr_ratio
            )  # water vapour content []
        else:
            self.qvap = qvap

        self.qcond = qcond  # liquid water content []
        self.WVEL = WVEL  # vertical velocity [m/s]
        self.UVEL = UVEL  # horizontal eastwards velocity [m/s]
        self.VVEL = VVEL  # horizontal northwards velocity [m/s]

    def generate_winds(self, ndims, ntime, THERMODATA):
        return constant_winds(ndims, ntime, THERMODATA, self.WVEL, self.UVEL, self.VVEL)

    def generate_thermo(self, gbxbounds, ndims, ntime):
        # shape_cen = ngridboxes * ntime = no. data for var defined at gridbox centers
        shape_cen = int(ntime * np.prod(ndims))
        THERMODATA = {
            "PRESS": np.full(shape_cen, self.PRESS),
            "TEMP": np.full(shape_cen, self.TEMP),
            "qvap": np.full(shape_cen, self.qvap),
            "qcond": np.full(shape_cen, self.qcond),
        }

        THERMODATA = self.generate_winds(ndims, ntime, THERMODATA)

        return THERMODATA


class SimpleThermo2DFlowField:
    """create thermodynamics that's constant in time
    with (P,T,qc) uniform throughout the domain with relative humidity
    = 0.95 below Zbase and a 2D (z,x) dependent flow field"""

    def __init__(
        self,
        config_filename,
        constants_filename,
        PRESS,
        TEMP,
        qvapmethod,
        qvapparams,
        Zbase,
        qcond,
        WMAX,
        Zlength,
        Xlength,
        VVEL,
    ):
        inputs = thermoinputsdict(config_filename, constants_filename)

        self.PRESS = PRESS  # pressure [Pa]
        self.TEMP = TEMP  # temperature [T]
        self.qcond = qcond  # liquid water content []

        # determine qvap [below, above] z (cloud) base
        self.Zbase = Zbase
        qvaps = qparams_to_qvap(qvapmethod, qvapparams, inputs["Mr_ratio"], PRESS, TEMP)
        self.qvap_below, self.qvap_above = qvaps

        self.WMAX = WMAX  # max velocities constant
        self.Zlength = Zlength  # wavelength of velocity modulation in z direction [m]
        self.Xlength = Xlength  # wavelength of velocity modulation in x direction [m]
        self.VVEL = VVEL  # horizontal (y) velocity

        self.RGAS_DRY = inputs["RGAS_DRY"]
        self.RGAS_V = inputs["RGAS_V"]
        self.RHO0 = inputs["RHO0"]

    def generate_qvap_profile(self, zfulls):
        qvap = np.where(zfulls >= self.Zbase, self.qvap_above, self.qvap_below)

        return qvap

    def rhotilda(self, ZCOORDS):
        """returns dimensionless rho_dry profile for use in stream function"""

        PRESS, TEMP = self.hydrostatic_adiabatic_thermo(ZCOORDS)

        RHO_DRY = PRESS / ((self.RGAS_DRY + self.qvap * self.RGAS_V) * TEMP)

        rhotilda = RHO_DRY / self.RHO0

        return rhotilda

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
        for VEL in ["WVEL", "UVEL", "VVEL"]:
            THERMODATA[VEL] = np.array([])

        if self.WMAX is not None:
            WVEL, UVEL = self.wvel_uvel_from_flowfield(THERMODATA, gbxbounds, ndims)
            THERMODATA["WVEL"] = np.tile(WVEL, ntime)
            THERMODATA["UVEL"] = np.tile(UVEL, ntime)

            if self.VVEL is not None:
                shape_yface = int((ndims[2] + 1) * ndims[0] * ndims[1] * ntime)
                THERMODATA["VVEL"] = np.full(shape_yface, self.VVEL)

        return THERMODATA

    def generate_thermo(self, gbxbounds, ndims, ntime):
        zfulls, xfulls, yfulls = rgrid.fullcoords_forallgridboxes(gbxbounds, ndims)

        qvap = self.generate_qvap_profile(zfulls)

        shape_cen = int(ntime * np.prod(ndims))
        THERMODATA = {
            "PRESS": np.full(shape_cen, self.PRESS),
            "TEMP": np.full(shape_cen, self.TEMP),
            "qvap": np.tile(qvap, ntime),
            "qcond": np.full(shape_cen, self.qcond),
        }

        THERMODATA = self.generate_winds(gbxbounds, ndims, ntime, THERMODATA)

        return THERMODATA


class ConstDryHydrostaticAdiabat:
    """create thermodynamics that's constant in time
    and in hydrostatic equillibrium with a dry adiabat
    accounting for the mass of water vapour in the air.
    Equations derived from Arabas et al. 2015 (sect 2.1)"""

    def __init__(
        self,
        config_filename,
        constants_filename,
        PRESSz0,
        THETA,
        qvapmethod,
        qvapparams,
        Zbase,
        qcond,
        WMAX,
        Zlength,
        Xlength,
        VVEL,
        moistlayer,
    ):
        inputs = thermoinputsdict(config_filename, constants_filename)

        ### parameters of profile ###
        self.PRESSz0 = PRESSz0  # pressure at z=0m [Pa]
        self.THETA = THETA  # (constant) dry potential temperature [K]
        self.qcond = qcond  # liquid mass mixing ratio []
        self.WMAX = WMAX  # max velocities constant
        self.Zlength = Zlength  # wavelength of velocity modulation in z direction [m]
        self.Xlength = Xlength  # wavelength of velocity modulation in x direction [m]
        self.VVEL = VVEL  # horizontal (y) velocity

        # determine qvap [below, above] z (cloud) base
        self.Zbase = Zbase
        self.qvapmethod, self.qvapparams = qvapmethod, qvapparams
        self.qvapz0 = qparams_to_qvap(
            qvapmethod, qvapparams, inputs["Mr_ratio"], self.PRESSz0, self.THETA
        )[0]
        self.moistlayer = moistlayer

        ### constants ###
        self.GRAVG = inputs["G"]
        self.CP_DRY = inputs["CP_DRY"]
        self.RGAS_DRY = inputs["RGAS_DRY"]
        self.RGAS_V = inputs["RGAS_V"]
        self.RC_DRY = self.RGAS_DRY / self.CP_DRY
        self.RCONST = 1 + self.qvapz0 * self.RGAS_V / self.RGAS_DRY
        self.P1000 = 100000  # P_1000 = 1000 hPa [Pa]
        self.CP0 = inputs["CP0"]
        self.RHO0 = inputs["RHO0"]
        self.Mr_ratio = inputs["Mr_ratio"]

        alpha = PRESSz0 / (self.RCONST * self.P1000)
        self.TEMPz0 = THETA * np.power(alpha, self.RC_DRY)  # temperature at z=0m [K]

        beta = (1 + self.qvapz0) / self.RCONST / self.RGAS_DRY
        self.RHOz0 = beta * self.PRESSz0 / self.TEMPz0

    def hydrostatic_adiabatic_profile(self, ZCOORDS):
        """returns *profile* of density (not the density itself!)
        rho = rhoprofile^((1-RC_DRY)/RC_DRY) = profile^pow"""

        pow = 1 / self.RC_DRY - 1

        Aa = (1 + self.qvapz0) * np.power(self.P1000, self.RC_DRY)
        Aa = self.THETA * self.RGAS_DRY / Aa
        Aconst = self.RCONST * np.power(Aa, (1 / (1 - self.RC_DRY)))

        RHOconst = -1 * self.GRAVG * self.RC_DRY / Aconst
        RHOprofile = np.power(self.RHOz0, 1 / pow) + RHOconst * ZCOORDS  # RHO^pow

        return RHOprofile, Aconst

    def hydrostatic_adiabatic_thermo(self, ZCOORDS):
        RHOprof, Aconst = self.hydrostatic_adiabatic_profile(ZCOORDS)

        # RHO = np.power(RHOprof, (1 / self.RC_DRY - 1))

        PRESS = Aconst * np.power(RHOprof, 1 / self.RC_DRY)

        TEMPconst = np.power(Aconst / (self.RCONST * self.P1000), self.RC_DRY)
        TEMP = self.THETA * TEMPconst * RHOprof

        return PRESS, TEMP

    def rhotilda(self, ZCOORDS):
        """returns dimensionless rho_dry profile for use in stream function"""

        PRESS, TEMP = self.hydrostatic_adiabatic_thermo(ZCOORDS)

        RHO_DRY = PRESS / ((self.RGAS_DRY + self.qvapz0 * self.RGAS_V) * TEMP)

        rhotilda = RHO_DRY / self.RHO0

        return rhotilda

    def wvel_uvel_from_flowfield(self, gbxbounds, ndims):
        zfaces = rgrid.coords_forgridboxfaces(gbxbounds, ndims, "z")[0]
        xfaces = rgrid.coords_forgridboxfaces(gbxbounds, ndims, "x")[1]
        rhotilda_zfaces = self.rhotilda(zfaces)
        rhotilda_xfaces = self.rhotilda(xfaces)

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
        for VEL in ["WVEL", "UVEL", "VVEL"]:
            THERMODATA[VEL] = np.array([])

        if self.WMAX is not None:
            WVEL, UVEL = self.wvel_uvel_from_flowfield(gbxbounds, ndims)
            THERMODATA["WVEL"] = np.tile(WVEL, ntime)
            THERMODATA["UVEL"] = np.tile(UVEL, ntime)

            if self.VVEL is not None:
                shape_yface = int((ndims[2] + 1) * ndims[0] * ndims[1] * ntime)
                THERMODATA["VVEL"] = np.full(shape_yface, self.VVEL)

        return THERMODATA

    def generate_qvap(self, zfulls, xfulls, PRESS, TEMP):
        qvaps = qparams_to_qvap(
            self.qvapmethod, self.qvapparams, self.Mr_ratio, PRESS, TEMP
        )
        qvap = np.where(zfulls < self.Zbase, qvaps[0], qvaps[1])

        if self.moistlayer:
            z1, z2 = self.moistlayer["z1"], self.moistlayer["z2"]
            x1, x2 = self.moistlayer["x1"], self.moistlayer["x2"]
            mlqvap = sratio2qvap(
                self.moistlayer["mlsratio"], PRESS, TEMP, self.Mr_ratio
            )
            moistregion = (
                (zfulls >= z1) & (zfulls < z2) & (xfulls >= x1) & (xfulls < x2)
            )
            qvap = np.where(moistregion, mlqvap, qvap)

        return qvap

    def generate_thermo(self, gbxbounds, ndims, ntime):
        zfulls, xfulls, yfulls = rgrid.fullcoords_forallgridboxes(gbxbounds, ndims)
        PRESS, TEMP = self.hydrostatic_adiabatic_thermo(zfulls)

        qvap = self.generate_qvap(zfulls, xfulls, PRESS, TEMP)

        shape_cen = int(ntime * np.prod(ndims))
        THERMODATA = {
            "PRESS": np.tile(PRESS, ntime),
            "TEMP": np.tile(TEMP, ntime),
            "qvap": np.tile(qvap, ntime),
            "qcond": np.full(shape_cen, self.qcond),
        }

        THERMODATA = self.generate_winds(gbxbounds, ndims, ntime, THERMODATA)

        return THERMODATA


class ConstHydrostaticLapseRates:
    """create thermodynamics that's constant in time
    and in hydrostatic equillibrium and following adiabats
    with constant lapse rates above/below zbase."""

    def __init__(
        self,
        config_filename,
        constants_filename,
        PRESS0,
        TEMP0,
        qvap0,
        Zbase,
        TEMPlapses,
        qvaplapses,
        qcond,
        WMAX,
        UVEL,
        VVEL,
        Wlength,
    ):
        self.PRESS0 = PRESS0  # surface pressure [Pa]
        self.TEMP0 = TEMP0  # surface temperature [T]
        self.qvap0 = qvap0  # surface water vapour content [Kg/Kg]
        self.Zbase = Zbase  # cloud base height [m]
        self.TEMPlapses = TEMPlapses  # temp lapse rates [below, above] Zbase [K km^-1]
        self.qvaplapses = (
            qvaplapses  # qvap lapse rates [below, above] Zbase [g/Kg km^-1]
        )

        self.qcond = qcond  # liquid water content [Kg/Kg]
        self.WMAX = WMAX  # vertical (coord3) velocity [m/s]
        self.UVEL = UVEL  # horizontal eastwards (coord1) velocity [m/s]
        self.VVEL = VVEL  # horizontal northwards (coord2) velocity [m/s]
        self.Wlength = Wlength  # [m] use constant W (Wlength=0.0), or sinusoidal 1-D profile below cloud base

        inputs = thermoinputsdict(config_filename, constants_filename)
        self.GRAVG = inputs["G"]
        self.RGAS_DRY = inputs["RGAS_DRY"]
        self.Mr_ratio = inputs["Mr_ratio"]

    def temp1(self, z):
        """note unit conversion of input lapse rates:
        templapse rate = -dT/dz [K km^-1]  -->  [K m^-1]"""

        temp1 = self.TEMP0 - self.TEMPlapses[0] / 1000 * z
        if np.any((temp1 <= 0.0)):
            raise ValueError("TEMP > 0.0K")
        return temp1

    def temp2(self, z):
        """note unit conversion of input lapse rates:
        templapse rate = -dT/dz [K km^-1]  -->  [K m^-1]"""

        T_Zbase = self.temp1(self.Zbase)  # TEMP at Zbase
        temp2 = T_Zbase - self.TEMPlapses[1] / 1000 * (z - self.Zbase)
        if np.any((temp2 <= 0.0)):
            raise ValueError("TEMP > 0.0K")
        return temp2

    def hydrostatic_pressure(self, P0, integral):
        exponent = -self.GRAVG / self.RGAS_DRY * integral
        return P0 * np.exp(exponent)

    def press1(self, z):
        """hydrostatic pressure for value z where z <= self.Zbase"""
        P0 = self.PRESS0
        integral = integrate.quad(lambda x: 1 / self.temp1(x), 0.0, z)[0]
        return self.hydrostatic_pressure(P0, integral)

    def press2(self, z):
        """hydrostatic pressure for value z where z > self.Zbase"""
        P0 = self.press1(self.Zbase)
        integral = integrate.quad(lambda x: 1 / self.temp2(x), self.Zbase, z)[0]
        return self.hydrostatic_pressure(P0, integral)

    def qvap1(self, z):
        """note unit conversion of input lapse rates:
        qvaplapse rate = -dqvap/dz [g/Kg km^-1]  -->  [m^-1]"""

        if self.qvaplapses[0] == "saturated":
            sratio = 1.001
            qvap1 = sratio2qvap(sratio, self.press2(z), self.temp2(z), self.Mr_ratio)
        else:
            qvap1 = self.qvap0 - self.qvaplapses[0] / 1e6 * z

        if np.any((qvap1 <= 0.0)):
            raise ValueError("TEMP > 0.0K")
        return qvap1

    def qvap2(self, z):
        """note unit conversion of input lapse rates:
        qvaplapse rate = -dqvap/dz [g/Kg km^-1]  -->  [m^-1]"""
        if self.qvaplapses[1] == "saturated":
            sratio = 1.001
            qvap2 = sratio2qvap(sratio, self.press2(z), self.temp2(z), self.Mr_ratio)
        else:
            qvap_Zbase = self.qvap1(self.Zbase)  # qvap at Zbase
            qvap2 = qvap_Zbase - self.qvaplapses[1] / 1e6 * (z - self.Zbase)

        if np.any((qvap2 <= 0.0)):
            raise ValueError("TEMP > 0.0K")
        return qvap2

    def below_above_zbase(self, zfulls, func1, func2):
        return np.where(zfulls <= self.Zbase, func1(zfulls), func2(zfulls))

    def below_above_zbase_qvap(self, zfulls):
        qvap = []
        for z in zfulls:
            if z < self.Zbase:
                qvap.append(self.qvap1(z))
            else:
                qvap.append(self.qvap2(z))

        return np.asarray(qvap)

    def below_above_zbase_pressure(self, zfulls):
        PRESS = []
        for z in zfulls:
            if z < self.Zbase:
                PRESS.append(self.press1(z))
            else:
                PRESS.append(self.press2(z))

        return np.asarray(PRESS)

    def hydrostatic_lapserates_thermo(self, zfulls):
        TEMP = self.below_above_zbase(zfulls, self.temp1, self.temp2)
        PRESS = self.below_above_zbase_pressure(zfulls)
        qvap = self.below_above_zbase_qvap(zfulls)

        return TEMP, PRESS, qvap

    def wvel_profile(self, gbxbounds, ndims, ntime):
        """returns updraught (w always >=0.0) sinusoidal
        profile with amplitude WMAX and wavelength 2*Wlength"""

        zfaces = rgrid.coords_forgridboxfaces(gbxbounds, ndims, "z")[0]
        WVEL = self.WMAX * np.sin(np.pi * zfaces / (2 * self.Wlength))

        WVEL[WVEL < 0.0] = 0.0

        return np.tile(WVEL, ntime)

    def generate_winds(self, gbxbounds, ndims, ntime, THERMODATA):
        THERMODATA = constant_winds(
            ndims, ntime, THERMODATA, self.WMAX, self.UVEL, self.VVEL
        )
        if self.Wlength > 0.0:
            THERMODATA["WVEL"] = self.wvel_profile(gbxbounds, ndims, ntime)

        return THERMODATA

    def generate_thermo(self, gbxbounds, ndims, ntime):
        zfulls, xfulls, yfulls = rgrid.fullcoords_forallgridboxes(gbxbounds, ndims)

        TEMP, PRESS, qvap = self.hydrostatic_lapserates_thermo(zfulls)

        shape_cen = int(ntime * np.prod(ndims))
        THERMODATA = {
            "PRESS": np.tile(PRESS, ntime),
            "TEMP": np.tile(TEMP, ntime),
            "qvap": np.tile(qvap, ntime),
            "qcond": np.full(shape_cen, self.qcond),
        }

        THERMODATA = self.generate_winds(gbxbounds, ndims, ntime, THERMODATA)

        return THERMODATA
