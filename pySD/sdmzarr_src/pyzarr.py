'''
----- CLEO -----
File: pyzarr.py
Project: sdmzarr_src
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

import thermodata
import pysetuptxt

def get_rawdataset(dataset):

  print("dataset: ", dataset)
  return xr.open_dataset(dataset, engine="zarr", consolidated=False)

def get_config(setuptxt):
  '''returns dictionary of configuration parameters
  read from from setup.txt file '''

  return pysetuptxt.config_dict(setuptxt)

def get_conts(setuptxt):
  '''returns dictionary of constants
  read from from setup.txt file '''

  return pysetuptxt.consts_dict(setuptxt)

def get_thermodata(dataset, ntime, ndims, consts):
  ''' returns a thermodynamic data in a dictionary. The value under 
  each key is the thermodynamics data in a 2D array 
  with dimensions [time, gridbox]. E.g. thermo["qvap"][:,0] gives the 
  timeseries of qvap for the 0th gridbox. thermo["qvap][0] gives 
  the qvap of all gridboxes at the 0th output time '''

  return thermodata.Thermodata(dataset, ntime, ndims, consts)
