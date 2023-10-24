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
from . import sddata

def get_rawdataset(dataset):

  print("dataset: ", dataset)
  return xr.open_dataset(dataset, engine="zarr", consolidated=False)

def raggedvar_fromzarr(ds, raggedcount, var):
  ''' returns ragged ak.Array dims [time, ragged]
  for a variable "var" in zarr ds '''

  return ak.unflatten(ds[var].values, raggedcount)

def var4d_fromzarr(ds, ntime, ndims, key):
  '''' returns 4D variable with dims
  [time, y, x, z] from zarr dataset "ds" '''
  
  reshape = [ntime] + list(ndims)

  return np.reshape(ds[key].values, reshape) 

def var3d_fromzarr(ds, ndims, key):
  '''' returns 3D variable with dims
  [y, x, z] from zarr dataset "ds" '''
  
  return np.reshape(ds[key].values, ndims) 

def get_thermodata(dataset, ntime, ndims, consts):
  ''' returns a thermodynamic data in a dictionary. The value under 
  each key is the thermodynamics data in a 2D array 
  with dimensions [time, gridbox]. E.g. thermo["qvap"][:,0] gives the 
  timeseries of qvap for the 0th gridbox. thermo["qvap][0] gives 
  the qvap of all gridboxes at the 0th output time '''

  return thermodata.Thermodata(dataset, ntime, ndims, consts)

def get_sddata(dataset):
  
  return sddata.SdData(dataset)

def get_gbxindex(dataset, ndims):

  if type(dataset) == str:
    dataset = get_rawdataset(dataset) 

  return var3d_fromzarr(dataset, ndims, "gbxindex")

def get_nsupers(dataset, ntime, ndims):

  if type(dataset) == str:
    dataset = get_rawdataset(dataset) 

  return var4d_fromzarr(dataset, ntime, ndims, "nsupers")