'''
----- CLEO -----
File: thermodata.py
Project: sdmout_src
Created Date: Tuesday 24th October 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
Last Modified: Monday 25th March 2024
Modified By: CB
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
Copyright (c) 2023 MPI-M, Clara Bayley
-----
File Description:
python class to handle thermodynamics and wind fields
data from SDM zarr store in cartesian
domain
'''

import numpy as np
import xarray as xr

from . import thermoeqns

def thermotryopen_dataset(dataset):

  if type(dataset) == str:
    print("thermodynamic fields dataset: ", dataset)
    return xr.open_dataset(dataset, engine="zarr", consolidated=False)
  else:
    return dataset

def thermovar4d_fromzarr(ds, reshape, key):
  '''' returns 4D variable with dims
  [time, y, x, z] from zarr dataset "ds" '''

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

    self.press_units = ds["press"].units # probably hecto pascals
    self.temp_units = ds["temp"].units # probably kelvin
    self.theta_units = ds["temp"].units # probably kelvin

  def potential_temp(self):
    ''' potential temperature, theta '''

    press = self.press*100 # convert from hPa to Pa
    theta = thermoeqns.dry_pot_temp(self.temp, press,
                                    self.qvap, self.consts)

    return theta

  def saturationpressure(self):
    ''' saturation pressure in hectoPascals '''

    psat = thermoeqns.saturation_pressure(self.temp) # [Pa]

    return  psat / 100 # [hPa]

  def vapourpressure(self):
    '''returns vapour pressure in hectoPascals '''

    p_pascals = self.press*100 # convert from hPa to Pa
    Mr_ratio = self.consts["Mr_ratio"]
    pvap = thermoeqns.vapour_pressure(p_pascals, self.qvap, Mr_ratio) # [Pa]

    return pvap / 100 # [hPa]

  def relative_humidity(self):
    ''' returns relative humidty and supersaturation '''

    p_pascals = self.press*100 # convert from hPa to Pa
    Mr_ratio = self.consts["Mr_ratio"]
    relh = thermoeqns.relative_humidity(p_pascals, self.temp,
                                        self.qvap, Mr_ratio)
    return relh

  def supersaturation(self):
    ''' returns relative humidty and supersaturation '''

    p_pascals = self.press*100 # convert from hPa to Pa
    Mr_ratio = self.consts["Mr_ratio"]
    supersat = thermoeqns.supersaturation(p_pascals, self.temp,
                                          self.qvap, Mr_ratio)
    return supersat

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
      err = "no known return provided for "+key+" key"
      raise ValueError(err)

class Winddata:

  def __init__(self, dataset, ntime, ndims, consts):

    ds = thermotryopen_dataset(dataset)

    self.consts = consts

    reshape = [ntime] + list(ndims)
    self.wvel = thermovar4d_fromzarr(ds, reshape, "wvel")
    self.uvel = thermovar4d_fromzarr(ds, reshape, "uvel")
    self.vvel = thermovar4d_fromzarr(ds, reshape, "vvel")

    self.wvel_units = ds["wvel"].units # probably hecto pascals
    self.uvel_units = ds["uvel"].units # probably hecto pascals
    self.vvel_units = ds["vvel"].units # probably hecto pascals

  def __getitem__(self, key):
    if key == "wvel":
      return self.wvel
    elif key == "uvel":
      return self.uvel
    elif key == "vvel":
      return self.vvel
    else:
      err = "no known return provided for "+key+" key"
      raise ValueError(err)
