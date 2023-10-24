'''
----- CLEO -----
File: pyzarr.py
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
functions to return zarr data in useful
formats for plotting e.g. for handling
ragged arrays of superdroplet attributes
'''

import numpy as np
import xarray as xr
import awkward as ak

from . import thermodata
from . import supersdata
from . import massmoms
from . import timedata

def get_rawdataset(dataset):
  return xr.open_dataset(dataset, engine="zarr", consolidated=False)

def raggedvar_fromzarr(ds, raggedcount, var):
  ''' returns ragged ak.Array dims [time, ragged]
  for a variable "var" in zarr ds '''
  if type(ds) == str:
    ds = get_rawdataset(ds) 

  return ak.unflatten(ds[var].values, raggedcount)

def var4d_fromzarr(ds, ntime, ndims, key):
  '''' returns 4D variable with dims
  [time, y, x, z] from zarr dataset "ds" '''
  if type(ds) == str:
    ds = get_rawdataset(ds) 

  reshape = [ntime] + list(ndims)
  return np.reshape(ds[key].values, reshape) 

def var3d_fromzarr(ds, ndims, key):
  '''' returns 3D variable with dims
  [y, x, z] from zarr dataset "ds" '''
  if type(ds) == str:
    ds = get_rawdataset(ds) 

  return np.reshape(ds[key].values, ndims) 

def get_thermodata(dataset, ntime, ndims, consts):
  ''' returns a thermodynamic data in a dictionary. The value under 
  each key is the thermodynamics data in a 2D array 
  with dimensions [time, gridbox]. E.g. thermo["qvap"][:,0] gives the 
  timeseries of qvap for the 0th gridbox. thermo["qvap][0] gives 
  the qvap of all gridboxes at the 0th output time '''

  return thermodata.Thermodata(dataset, ntime, ndims, consts)

def get_supers(dataset):
  return supersdata.SupersData(dataset)

def get_rainsupers(dataset):
  return supersdata.RainSupers(dataset)

def get_time(dataset):
  return timedata.Time(dataset)

def get_massmoms(dataset, ntime, ndims):
  return massmoms.MassMoms(dataset, ntime, ndims)

def get_gbxindex(dataset, ndims):
  return var3d_fromzarr(dataset, ndims, "gbxindex")

def get_totnsupers(dataset):
  
  if type(dataset) == str:
    dataset = get_rawdataset(dataset) 
  try:
    return dataset["rgd_totnsupers"]
  except:
    return dataset["totnsupers"]

def get_nsupers(dataset, ntime, ndims):
  return var4d_fromzarr(dataset, ntime, ndims, "nsupers")

def get_nrainsupers(dataset, ntime, ndims):
  return var4d_fromzarr(dataset, ntime, ndims, "nrainsupers")