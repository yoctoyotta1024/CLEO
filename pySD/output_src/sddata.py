'''
----- CLEO -----
File: sddata.py
Project: output_src
Created Date: Tuesday 24th October 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
Last Modified: Tuesday 24th October 2023
Modified By: CB
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
Copyright (c) 2023 MPI-M, Clara Bayley
-----
File Description:
'''

import numpy as np
import xarray as xr
import awkward as ak

class SdData:
  
  def __init__(self, dataset):
    
    ds = self.tryopen_dataset(dataset) 
    rgdcount = ds["rgd_totnsupers"].values # ragged count variable
     
    self.sdindex = self.tryvar(ds, rgdcount, "sdId")
    self.sd_gbxindex = self.tryvar(ds, rgdcount, "sdgbxindex")
    self.xi = self.tryvar(ds, rgdcount, "xi")
    self.radius = self.tryvar(ds, rgdcount, "radius")
    self.m_sol = self.tryvar(ds, rgdcount, "msol")

    self.coord3 = self.tryvar(ds, rgdcount, "coord3")
    self.coord1 = self.tryvar(ds, rgdcount, "coord1")
    self.coord2 = self.tryvar(ds, rgdcount, "coord2")

    self.radius_units = self.tryunits(ds, "radius") # probably microns ie. 'micro m'
    self.m_sol_units = self.tryunits(ds, "msol") # probably gramms
    self.coord3_units = self.tryunits(ds, "coord3") # probably meters
    self.coord1_units = self.tryunits(ds, "coord1") # probably meters
    self.coord2_units = self.tryunits(ds, "coord2") # probably meters
  
  def tryopen_dataset(dataset):
    
    if type(dataset) == str:
      return xr.open_dataset(dataset, engine="zarr", consolidated=False) 
    else:
      return dataset
    
  def tryvar(ds, raggedcount, var):
    ''' attempts to return variable in form of ragged array
    (ak.Array) with dims [time, raggedcount]
    for a variable "var" in xarray dataset 'ds'.
    If attempt fails, returns empty array instead '''
    try:
      return ak.unflatten(ds[var].values, raggedcount)
    except:
      return ak.Array([])

  def tryunits(ds, var):
    ''' attempts to return the units of a variable 
    in xarray dataset 'ds'. If attempt fails, returns null '''
    try:
      return ds[var].units
    except:
      return ""

  def __getitem__(self, key):

    if key == "sdId":
      return self.sdId
    elif key == "sdgbxindex":
      return self.sd_gbxindex
    elif key == "xi":
      return self.xi
    elif key == "radius":
      return self.radius
    elif key == "msol":
      return self.m_sol
    elif key == "coord3":
      return self.coord3
    elif key == "coord1":
      return self.coord1
    elif key == "coord2":
      return self.coord2
    else:
      err = "no known return provided for "+key+" key"
      raise ValueError(err)
