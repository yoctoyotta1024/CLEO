"""
Copyright (c) 2024 MPI-M, Clara Bayley


----- CLEO -----
File: thermodata.py
Project: sdmout_src
Created Date: Tuesday 24th October 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
python class to handle thermodynamics and wind fields data from SDM zarr store in cartesian domain
"""

import numpy as np
import xarray as xr
from pathlib import Path

from . import thermoeqns


def thermotryopen_dataset(dataset):
    if isinstance(dataset, str) or isinstance(dataset, Path):
        print("thermodynamic fields dataset: ", dataset)
        return xr.open_dataset(dataset, engine="zarr", consolidated=False)
    else:
        return dataset


def thermovar4d_fromzarr(ds, reshape, key):
    """' returns 4D variable with dims
    [time, y, x, z] from zarr dataset "ds" """

    return np.reshape(ds[key].values, reshape)


class Thermodata:
    def __init__(self, dataset, ntime, ndims, consts):
        ds = thermotryopen_dataset(dataset)

        self.consts = consts

        reshape = [ntime] + list(ndims)
        self.press = thermovar4d_fromzarr(ds, reshape, "press")
        self.temp = thermovar4d_fromzarr(ds, reshape, "temp")
        self.qvap = thermovar4d_fromzarr(ds, reshape, "qvap")
        self.qcond = thermovar4d_fromzarr(ds, reshape, "qcond")
        self.theta = self.potential_temp()

        self.press_units = ds["press"].units  # assumed probably hecto-pascals
        self.temp_units = ds["temp"].units  # assumed probably kelvin
        self.qvap_units = ds["qvap"].units  # assumed probably g/Kg
        self.qcond_units = ds["qcond"].units  # assumed probably g/Kg
        self.theta_units = ds["temp"].units  # assumed probably kelvin

    def potential_temp(self):
        """potential temperature, theta"""

        press = self.press * 100  # convert from hPa to Pa
        qvap = self.qvap / 1000  # convert g/Kg to Kg/Kg

        return thermoeqns.dry_pot_temp(self.temp, press, qvap, self.consts)

    def saturationpressure(self):
        """saturation pressure in hectoPascals"""

        psat = thermoeqns.saturation_pressure(self.temp)  # [Pa]

        return psat / 100  # [hPa]

    def vapourpressure(self):
        """returns vapour pressure in hectoPascals"""

        p_pascals = self.press * 100  # convert from hPa to Pa
        qvap = self.qvap / 1000  # convert g/Kg to Kg/Kg
        Mr_ratio = self.consts["Mr_ratio"]
        pvap = thermoeqns.vapour_pressure(p_pascals, qvap, Mr_ratio)  # [Pa]

        return pvap / 100  # [hPa]

    def relative_humidity(self):
        """returns relative humidty and supersaturation"""

        p_pascals = self.press * 100  # convert from hPa to Pa
        qvap = self.qvap / 1000  # convert g/Kg to Kg/Kg
        Mr_ratio = self.consts["Mr_ratio"]

        return thermoeqns.relative_humidity(p_pascals, self.temp, qvap, Mr_ratio)

    def supersaturation(self):
        """returns relative humidty and supersaturation"""

        p_pascals = self.press * 100  # convert from hPa to Pa
        qvap = self.qvap / 1000  # convert g/Kg to Kg/Kg
        Mr_ratio = self.consts["Mr_ratio"]
        return thermoeqns.supersaturation(p_pascals, self.temp, qvap, Mr_ratio)

    def __getitem__(self, key):
        if key == "press":
            return self.press
        elif key == "temp":
            return self.temp
        elif key == "qvap":
            return self.qvap
        elif key == "qcond":
            return self.qcond
        elif key == "theta":
            return self.theta
        elif key == "relh":
            return self.relative_humidity()
        elif key == "supersat":
            return self.supersaturation()
        else:
            err = "no known return provided for " + key + " key"
            raise ValueError(err)


class Winddata:
    def __init__(self, dataset, ntime, ndims, consts):
        ds = thermotryopen_dataset(dataset)

        self.consts = consts

        reshape = [ntime] + list(ndims)
        self.wvel = thermovar4d_fromzarr(ds, reshape, "wvel")
        self.uvel = thermovar4d_fromzarr(ds, reshape, "uvel")
        self.vvel = thermovar4d_fromzarr(ds, reshape, "vvel")

        self.wvel_units = ds["wvel"].units  # probably hecto pascals
        self.uvel_units = ds["uvel"].units  # probably hecto pascals
        self.vvel_units = ds["vvel"].units  # probably hecto pascals

    def __getitem__(self, key):
        if key == "wvel":
            return self.wvel
        elif key == "uvel":
            return self.uvel
        elif key == "vvel":
            return self.vvel
        else:
            err = "no known return provided for " + key + " key"
            raise ValueError(err)
