"""
Copyright (c) 2025 MPI-M, Clara Bayley


----- CLEO -----
File: supersdata.py
Project: sdmout_src
Created Date: Tuesday 24th October 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
!!!NOTE!!! DEPRECIATED SINCE V0.40.0, REPLACED BY SUPERDROPS MODULE !!!!!!
python class to handle superdroplet attributes data from SDM zarr store in ragged array format
"""


import numpy as np
import xarray as xr
import awkward as ak
import warnings
from pathlib import Path
from typing import Optional

warnings.warn(
    "supersdata module is depreciated since CLEO v0.40.0,\nreplaced by superdrops"
)


class SuperdropProperties:
    """Contains attributes common to all superdroplets and functions
    for calculating derived ones"""

    def __init__(self, consts):
        """Common attributes shared by superdroplets"""

        # density of liquid in droplets (=density of water at 300K) [Kg/m^3]
        self.RHO_L = consts["RHO_L"]

        # density of (dry) solute [Kg/m^3]
        self.RHO_SOL = consts["RHO_SOL"]

        # Mr of solute [g/mol]
        self.MR_SOL = consts["MR_SOL"]

        # degree ionic dissociation (van't Hoff factor)
        self.IONIC = consts["IONIC"]

        self.print_properties()

    def print_properties(self):
        print("---- Superdrop Properties -----")
        print("RHO_L =", self.RHO_L, "Kg/m^3")
        print("RHO_SOL =", self.RHO_SOL, "Kg/m^3")
        print("MR_SOL =", self.MR_SOL, "Kg/mol")
        print("IONIC =", self.IONIC)
        print("-------------------------------")

    def rhoeff(self, r, msol):
        """calculates effective density [g m^-3] of
        droplet such that mass_droplet, m = 4/3*pi*r^3 * rhoeff
        taking into account mass of liquid and mass of
        solute assuming solute occupies volume it
        would given its (dry) density, RHO_SOL."""

        msol = msol / 1000  # convert from grams to Kg
        r = r / 1e6  # convert microns to m

        solfactor = 3 * msol / (4.0 * np.pi * (r**3))
        rhoeff = self.RHO_L + solfactor * (1 - self.RHO_L / self.RHO_SOL)

        return rhoeff * 1000  # [g/m^3]

    def vol(self, r):
        """volume of droplet [m^3]"""

        r = r / 1e6  # convert microns to m

        return 4.0 / 3.0 * np.pi * r**3

    def mass(self, r, msol):
        """
        total mass of droplet (water + (dry) areosol) [g],
        m =  4/3*pi*rho_l**3 + msol(1-rho_l/rho_sol)
        ie. m = 4/3*pi*rhoeff*R**3
        """

        msol = msol / 1000  # convert from grams to Kg
        r = r / 1e6  # convert microns to m

        msoleff = msol * (1 - self.RHO_L / self.RHO_SOL)  # effect of solute on mass
        m = msoleff + 4 / 3.0 * np.pi * (r**3) * self.RHO_L

        return m * 1000  # [g]

    def m_water(self, r, msol):
        """mass of only water in droplet [g]"""

        msol = msol / 1000  # convert msol from grams to Kg
        r = r / 1e6  # convert microns to m

        v_sol = msol / self.RHO_SOL
        v_w = 4 / 3.0 * np.pi * (r**3) - v_sol

        return self.RHO_L * v_w * 1000  # [g]


class SupersData(SuperdropProperties):
    def __init__(self, dataset: Path, consts: Optional[dict] = None):
        if consts is not None:
            SuperdropProperties.__init__(self, consts)
        else:
            warnings.warn("No superdroplet properties attached to superdroplet dataset")

        ds = self.tryopen_dataset(dataset)
        raggedcount = ds["raggedcount"]  # ragged count variable

        self.sdId = self.tryvar(ds, raggedcount, "sdId")
        self.sdgbxindex = self.tryvar(ds, raggedcount, "sdgbxindex")
        self.xi = self.tryvar(ds, raggedcount, "xi")
        self.radius = self.tryvar(ds, raggedcount, "radius")
        self.msol = self.tryvar(ds, raggedcount, "msol")

        self.coord3 = self.tryvar(ds, raggedcount, "coord3")
        self.coord1 = self.tryvar(ds, raggedcount, "coord1")
        self.coord2 = self.tryvar(ds, raggedcount, "coord2")

        self.radius_units = self.tryunits(ds, "radius")  # probably microns
        self.msol_units = self.tryunits(ds, "msol")  # probably gramms
        self.coord3_units = self.tryunits(ds, "coord3")  # probably meters
        self.coord1_units = self.tryunits(ds, "coord1")  # probably meters
        self.coord2_units = self.tryunits(ds, "coord2")  # probably meters

    def tryopen_dataset(self, dataset: Path):
        if isinstance(dataset, str) or isinstance(dataset, Path):
            print("supers dataset: ", dataset)
            return xr.open_dataset(dataset, engine="zarr", consolidated=False)
        else:
            return dataset

    def tryvar(self, ds: xr.Dataset, raggedcount: xr.DataArray, var: str):
        """attempts to return variable in form of ragged array (ak.Array) with dims
        [time, raggedcount] for a variable "var" in xarray dataset 'ds'.
        If attempt fails, returns empty array instead"""
        try:
            return ak.unflatten(ds[var].values, raggedcount.values)
        except KeyError:
            return ak.Array([])

    def tryunits(self, ds: xr.Dataset, var: str):
        """attempts to return the units of a variable
        in xarray dataset 'ds'. If attempt fails, returns null"""
        try:
            return ds[var].units
        except KeyError:
            return ""

    def __getitem__(self, key):
        if key == "sdId":
            return self.sdId
        elif key == "sdgbxindex":
            return self.sdgbxindex
        elif key == "xi":
            return self.xi
        elif key == "radius":
            return self.radius
        elif key == "msol":
            return self.msol
        elif key == "coord3":
            return self.coord3
        elif key == "coord1":
            return self.coord1
        elif key == "coord2":
            return self.coord2
        else:
            err = "no known return provided for " + key + " key"
            raise ValueError(err)


class RainSupers(SuperdropProperties):
    def __init__(self, sddata, consts, rlim=40):
        """return data for (rain)drops with radii > rlim.
        Default minimum raindrops size is rlim=40microns"""

        if not isinstance(sddata, SupersData):
            sddata = SupersData(dataset=sddata, consts=consts)

        israin = sddata.radius >= rlim  # ak array True for raindrops

        self.totnsupers_rain = ak.num(israin[israin is True])
        self.sdId = sddata.sdId[israin]
        self.sdgbxindex = sddata.sdgbxindex[israin]
        self.xi = sddata.xi[israin]
        self.radius = sddata.radius[israin]
        self.msol = sddata.msol[israin]

        if np.any(sddata.coord3):
            self.coord3 = sddata.coord3[israin]
            if np.any(sddata.coord1):
                self.coord1 = sddata.coord1[israin]
                if np.any(sddata.coord2):
                    self.coord2 = sddata.coord2[israin]

    def __getitem__(self, key):
        if key == "totnsupers_rain":
            return self.totnsupers_rain
        elif key == "sdId":
            return self.sdId
        elif key == "sdgbxindex":
            return self.sdgbxindex
        elif key == "xi":
            return self.xi
        elif key == "radius":
            return self.radius
        elif key == "msol":
            return self.msol
        elif key == "coord3":
            return self.coord3
        elif key == "coord1":
            return self.coord1
        elif key == "coord2":
            return self.coord2
        else:
            err = "no known return provided for " + key + " key"
            raise ValueError(err)
