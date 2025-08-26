"""
----- CLEO -----
File: thermoeqns.py
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
equations in python for calculating thermodynamic
variables e.g. potential temperature
"""

import numpy as np


def vapour_pressure(press, qv, Mr_ratio):
    pv = qv * press / (Mr_ratio + qv)

    return pv


def relative_humidity(p, temp, qv, Mr_ratio):
    pv = vapour_pressure(p, qv, Mr_ratio)
    psat = saturation_pressure(temp)
    relh = pv / psat

    return relh


def supersaturation(p, temp, qv, Mr_ratio):
    pv = vapour_pressure(p, qv, Mr_ratio)
    psat = saturation_pressure(temp)

    qsat = Mr_ratio * psat / (p - pv)
    supersat = qv / qsat - 1

    return supersat


def saturation_pressure(TEMP):
    """Calculate the equilibrium vapor pressure
    of water over liquid water ie. the
    saturation pressure (psat) given TEMP /K.
    Equation taken from Bjorn Steven's
    "make_tetens" python function from his module
    "moist_thermodynamics.saturation_vapour_pressures"
    available on gitlab. Original paper
    "Murray, F. W. On the Computation of Saturation
    Vapor Pressure. Journal of Applied Meteorology
    and Climatology 6, 203–204 (1967)." """

    A = 17.4146  # constants from Bjorn Gitlab originally from paper
    B = 33.639  # ditto
    TREF = 273.16  # Triple point temperature [K] of water
    PREF = 611.655  # Triple point pressure [Pa] of water

    if np.any(TEMP <= 0):
        err = "Temperature must be larger than 0K."
        raise ValueError(err)

    expothis = A * (TEMP - TREF) / (TEMP - B)

    return PREF * np.exp(expothis)  # dimensionless psat /Pa


def saturation_pressure_murphy_koop(TEMP):
    """Calculate the equilibrium vapor pressure
    of water over liquid water ie. the
    saturation pressure (psat [Pa]). Equation taken from
    typhon.physics.thermodynamics.e_eq_water_mk."""

    if np.any(TEMP <= 0):
        err = "Temperature must be larger than 0K."
        raise ValueError(err)

    lnpsat = (
        54.842763  # ln(psat) [Pa]
        - 6763.22 / TEMP
        - 4.21 * np.log(TEMP)
        + 0.000367 * TEMP
        + np.tanh(0.0415 * (TEMP - 218.8))
        * (53.878 - 1331.22 / TEMP - 9.44523 * np.log(TEMP) + 0.014025 * TEMP)
    )

    return np.exp(lnpsat)  # psat [Pa]


def dry_pot_temp(Temp, P, qv, consts):
    """calculate potential Temperature [K]
    assuming moist (unsaturated) air with
    vapour content qv"""

    Cpdry = consts["CP_DRY"]
    Cpv = consts["CP_V"]
    Rgasdry = consts["RGAS_DRY"]
    Rgasv = consts["RGAS_V"]

    Cp = Cpdry * (1 + qv * Cpv / Cpdry) / (1 + qv)
    Rgas = Rgasdry * (1 + qv * Rgasv / Rgasdry) / (1 + qv)

    Theta = Temp * (P[0] / P) ** (Rgas / Cp)

    return Theta


def dry_adiabat(p, temp, consts):
    gamma = consts["RGAS_DRY"] / consts["CP_DRY"]
    dry_adia = temp[0] * (p / p[0]) ** gamma  # dry adiabatic temp

    dry_adia_theta = dry_adia * (p[0] / p) ** gamma  # dry adiabatic theta (=const)

    return dry_adia, dry_adia_theta


def specific_moist_static_energy(Z, Temp, qv, consts):
    """calculate the specific mass moist static
    energy, mse, J/Kg. (ie. mse per unit mass)."""

    GRAVG = consts["G"]  # [m/s^2]
    Latent_v = consts["LATENT_V"]  # [J/Kg]
    Cp_dry = consts["CP_DRY"]  # [J/Kg/K]

    return GRAVG * Z + Latent_v * qv + Cp_dry * Temp
