'''
----- CLEO -----
File: thermodata.py
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
python class to handle thermodynamic
data from SDM zarr store in cartesian
domain
'''

import numpy as np
import xarray as xr

import thermoeqns

def getds(dataset):

  return xr.open_dataset(dataset, engine="zarr", consolidated=False)

