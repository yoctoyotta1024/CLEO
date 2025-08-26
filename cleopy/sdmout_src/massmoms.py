"""
Copyright (c) 2024 MPI-M, Clara Bayley


----- CLEO -----
File: massmoms.py
Project: sdmout_src
Created Date: Tuesday 24th October 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
python class to handle mass moments data from SDM zarr store in cartesian domain
"""

import numpy as np
import xarray as xr
from pathlib import Path


class MassMoms:
    def __init__(self, dataset, ntime, ndims, lab=""):
        ds = self.tryopen_dataset(dataset)
        reshape = [ntime] + list(ndims)

        try:
            self.nsupers = self.var4d_fromzarr(
                ds, reshape, "n" + lab + "supers"
            )  # number of superdroplets in gbxs over time
        except KeyError:
            self.nsupers = np.array([])

        self.mom0 = self.var4d_fromzarr(
            ds, reshape, "massmom0" + lab
        )  # number of droplets in gbxs over time
        self.mom1 = self.var4d_fromzarr(
            ds, reshape, "massmom1" + lab
        )  # total mass of droplets in gbxs over time
        self.mom2 = self.var4d_fromzarr(
            ds, reshape, "massmom2" + lab
        )  # 2nd mass moment of droplets (~reflectivity)
        self.effmass = self.effective_mass()

        self.mom1_units = ds["massmom1"].units  # probably grams
        self.mom2_units = ds["massmom2"].units  # probably grams^2
        self.effmass_units = (
            ds["massmom2"].units + "/" + ds["massmom1"].units
        )  # probably grams

    def tryopen_dataset(self, dataset):
        if isinstance(dataset, str) or isinstance(dataset, Path):
            print("mass moments dataset: ", dataset)
            return xr.open_dataset(dataset, engine="zarr", consolidated=False)
        else:
            return dataset

    def var4d_fromzarr(self, ds, reshape, key):
        """' returns 4D variable with dims
        [time, y, x, z] from zarr dataset "ds" """

        return np.reshape(ds[key].values, reshape)

    def effective_mass(self):
        """effective mass of droplets"""
        return self.mom2 / self.mom1

    def __getitem__(self, key):
        if key == "nsupers":
            return self.nsupers
        elif key == "mom0":
            return self.mom0
        elif key == "mom1":
            return self.mom1
        elif key == "mom2":
            return self.mom2
        elif key == "effmass":
            return self.effmass
        else:
            err = "no known return provided for " + key + " key"
            raise ValueError(err)
