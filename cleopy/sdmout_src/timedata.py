"""
Copyright (c) 2024 MPI-M, Clara Bayley


----- CLEO -----
File: timedata.py
Project: sdmout_src
Created Date: Tuesday 7th May 2024
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
python class to get time data from SDM zarr store and return in various formats
"""

import xarray as xr
from pathlib import Path


class Time:
    def __init__(self, dataset):
        timeds = self.tryopen_time(dataset)  # time from dataset

        if timeds.units != "s":
            raise ValueError("units of time in dataset must be seconds")

        self.secs = timeds.values
        self.mins = self.secs / 60
        self.hrs = self.secs / 60 / 60

    def tryopen_time(self, dataset):
        """returns time variable (with metadata) from xarray dataset"""

        if isinstance(dataset, str) or isinstance(dataset, Path):
            print("time from dataset: ", dataset)
            timeds = xr.open_dataset(dataset, engine="zarr", consolidated=False)["time"]
        else:
            timeds = dataset["time"]

        return timeds

    def __getitem__(self, key):
        if key == "secs":
            return self.secs
        elif key == "mins":
            return self.mins
        elif key == "hrs":
            return self.hrs
        else:
            err = "no known return provided for " + key + " key"
            raise ValueError(err)
