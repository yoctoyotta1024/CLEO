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

def get_rawdataset(dataset):

  print("dataset: ", dataset)
  return xr.open_dataset(dataset, engine="zarr", consolidated=False)
