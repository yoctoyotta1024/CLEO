"""
Copyright (c) 2024 MPI-M, Clara Bayley


----- CLEO -----
File: thermogen.py
Project: thermobinary_src
Created Date: Monday 16th October 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
Various ways of generating fields for Temp, Pressure, Qvap and Qcond to read into CLEO
"""

import numpy as np
from scipy import integrate
from .. import cxx2py
from .create_thermodynamics import thermoinputsdict
from .windsgen import DryAdiabat2DFlowField
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


class ConstUniformThermo:
    """create thermodynamics that's constant in time and uniform throughout the domain"""

    def __init__(
        self,
        PRESS,
        TEMP,
        qvap,
        qcond,
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

    def generate_thermo(self, gbxbounds, ndims, ntime):
        shape_cen = int(
            ntime * np.prod(ndims)
        )  # = no. data for var defined at gridbox centers

        THERMODATA = {
            "PRESS": np.full(shape_cen, self.PRESS),
            "TEMP": np.full(shape_cen, self.TEMP),
            "qvap": np.full(shape_cen, self.qvap),
            "qcond": np.full(shape_cen, self.qcond),
        }

        return THERMODATA


class Simple2TierRelativeHumidity:
    """create thermodynamics that's constant in time with (P,T,qc) uniform throughout the domain
    and with relative humidity uniform above and below Zbase"""

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
    ):
        inputs = thermoinputsdict(config_filename, constants_filename)

        self.PRESS = PRESS  # pressure [Pa]
        self.TEMP = TEMP  # temperature [T]
        self.qcond = qcond  # liquid water content []

        # determine qvap [below, above] z (cloud) base
        self.Zbase = Zbase
        qvaps = qparams_to_qvap(qvapmethod, qvapparams, inputs["Mr_ratio"], PRESS, TEMP)
        self.qvap_below, self.qvap_above = qvaps

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

    def generate_thermo(self, gbxbounds, ndims, ntime):
        zfulls, xfulls, yfulls = rgrid.fullcoords_forallgridboxes(gbxbounds, ndims)

        qvap = self.generate_qvap_profile(zfulls)

        shape_cen = int(
            ntime * np.prod(ndims)
        )  # = no. data for var defined at gridbox centers
        THERMODATA = {
            "PRESS": np.full(shape_cen, self.PRESS),
            "TEMP": np.full(shape_cen, self.TEMP),
            "qvap": np.tile(qvap, ntime),
            "qcond": np.full(shape_cen, self.qcond),
        }

        return THERMODATA


class DryHydrostaticAdiabatic2TierRelH:
    """create thermodynamics that's constant in time and in hydrostatic equillibrium with
    a dry adiabat accounting for the mass of water vapour in the air.
    Equations derived from Arabas et al. 2015 (sect 2.1).
    Relative humidity like for Simple2TierRelativeHumidity exceptional moist layer possible
    to add within in a certain height range"""

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
        moistlayer,
    ):
        inputs = thermoinputsdict(config_filename, constants_filename)

        ### parameters of profile ###
        self.PRESSz0 = PRESSz0  # pressure at z=0m [Pa]
        self.THETA = THETA  # (constant) dry potential temperature [K]
        self.qcond = qcond  # liquid mass mixing ratio []

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
        TEMPz0 = THETA * np.power(alpha, self.RC_DRY)  # temperature at z=0m [K]
        beta = (1 + self.qvapz0) / self.RCONST / self.RGAS_DRY
        self.RHOz0 = beta * self.PRESSz0 / TEMPz0

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

    def create_default_windsgen(self, WMAX, Zlength, Xlength, VVEL):
        return DryAdiabat2DFlowField(WMAX, Zlength, Xlength, VVEL, self)

    def generate_thermo(self, gbxbounds, ndims, ntime):
        zfulls, xfulls, yfulls = rgrid.fullcoords_forallgridboxes(gbxbounds, ndims)
        PRESS, TEMP = self.hydrostatic_adiabatic_thermo(zfulls)

        qvap = self.generate_qvap(zfulls, xfulls, PRESS, TEMP)

        shape_cen = int(
            ntime * np.prod(ndims)
        )  # = no. data for var defined at gridbox centers
        THERMODATA = {
            "PRESS": np.tile(PRESS, ntime),
            "TEMP": np.tile(TEMP, ntime),
            "qvap": np.tile(qvap, ntime),
            "qcond": np.full(shape_cen, self.qcond),
        }

        return THERMODATA


class HydrostaticLapseRates:
    """create thermodynamics that's constant in time and in hydrostatic equillibrium and
    following temperature and qvap (adiabats) with constant lapse rates above/below zbase.
    Qcond is uniform and constant."""

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

    def generate_thermo(self, gbxbounds, ndims, ntime):
        zfulls, xfulls, yfulls = rgrid.fullcoords_forallgridboxes(gbxbounds, ndims)

        TEMP, PRESS, qvap = self.hydrostatic_lapserates_thermo(zfulls)

        shape_cen = int(
            ntime * np.prod(ndims)
        )  # = no. data for var defined at gridbox centers
        THERMODATA = {
            "PRESS": np.tile(PRESS, ntime),
            "TEMP": np.tile(TEMP, ntime),
            "qvap": np.tile(qvap, ntime),
            "qcond": np.full(shape_cen, self.qcond),
        }

        return THERMODATA
