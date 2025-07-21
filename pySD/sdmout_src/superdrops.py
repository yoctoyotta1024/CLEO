"""
Copyright (c) 2025 MPI-M, Clara Bayley


----- CLEO -----
File: superdrops.py
Project: sdmout_src
Created Date: Thursday 29th May 2025
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
python class(es) to handle superdroplet attributes data from SDM zarr store in ragged array format.
Non-compatible update on supersdata and sdtracing modules
"""


import awkward as ak
import numpy as np
import random
import warnings
import xarray as xr

from pathlib import Path


class SuperdropProperties:
    """Contains attributes common to all superdroplets and functions
    for calculating derived ones"""

    def __init__(self, consts: dict, is_print=False):
        """Common attributes shared by superdroplets:
        RHO_L = density of liquid in droplets (=density of water at 300K) [Kg/m^3]
        RHO_SOL = density of (dry) solute [Kg/m^3]
        MR_SOL = Mr of solute [g/mol]
        IONIC = degree ionic dissociation (van't Hoff factor)
        """
        # properties: value
        self._props = {
            "RHO_L": consts["RHO_L"],
            "RHO_SOL": consts["RHO_SOL"],
            "MR_SOL": consts["MR_SOL"],
            "IONIC": consts["IONIC"],
            "RHO_L_units": "Kg m^{-3}",
            "RHO_SOL_units": "Kg m^{-3}",
            "MR_SOL_units": "g mol^{-1}",
        }

        if is_print:
            self.print_properties()

    def get_props(self, var: str, units: bool):
        if not units:
            return self._props[var]
        else:
            return self._props[var], self._props[var + "_units"]

    def RHO_L(self, units=False):
        return self.get_props("RHO_L", units=units)

    def RHO_SOL(self, units=False):
        return self.get_props("RHO_SOL", units=units)

    def MR_SOL(self, units=False):
        return self.get_props("MR_SOL", units=units)

    def IONIC(self, units=False):
        if units:
            print("degree ionic dissociation (van't Hoff factor) has no units")
        return self.get_props("IONIC", units=False)

    def __getitem__(self, key):
        return self._props[key]

    def print_properties(self):
        print("---- Superdrop Properties -----")
        print("RHO_L =", self.RHO_L())
        print("RHO_SOL =", self.RHO_SOL())
        print("MR_SOL =", self.MR_SOL())
        print("IONIC =", self.IONIC())
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
        rhoeff = self._props["RHO_L"] + solfactor * (
            1 - self._props["RHO_L"] / self._props["RHO_SOL"]
        )

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

        msoleff = msol * (
            1 - self._props["RHO_L"] / self._props["RHO_SOL"]
        )  # effect of solute on mass
        m = msoleff + 4 / 3.0 * np.pi * (r**3) * self._props["RHO_L"]

        return m * 1000  # [g]

    def m_water(self, r, msol):
        """mass of only water in droplet [g]"""

        msol = msol / 1000  # convert msol from grams to Kg
        r = r / 1e6  # convert microns to m

        v_sol = msol / self._props["RHO_SOL"]
        v_w = 4 / 3.0 * np.pi * (r**3) - v_sol

        return self._props["RHO_L"] * v_w * 1000  # [g]


class Superdrops(SuperdropProperties):
    def __init__(self, dataset, consts=None, is_print=False):
        if isinstance(dataset, Superdrops):
            SuperdropProperties.__init__(self, dataset._props, is_print=is_print)
            self._variables = dataset._variables
            self._raw_data = dataset._raw_data
            self._units = dataset._units

        elif (
            isinstance(dataset, xr.Dataset)
            or isinstance(dataset, str)
            or isinstance(dataset, Path)
        ):
            if consts is not None:
                SuperdropProperties.__init__(self, consts, is_print=is_print)
            else:
                warnings.warn(
                    "No superdroplet properties attached to superdroplet dataset"
                )

            ds = self.tryopen_dataset(dataset)
            raggedcount = ds["raggedcount"]  # ragged count variable
            self._variables = (
                "sdId",
                "sdgbxindex",
                "xi",
                "radius",
                "msol",
                "coord3",
                "coord1",
                "coord2",
            )
            self._raw_data = {
                var: self.tryvar(ds, raggedcount, var) for var in self._variables
            }
            self._units = {
                "radius_units": self.tryunits(ds, "radius"),  # probably microns
                "msol_units": self.tryunits(ds, "msol"),  # probably gramms
                "coord3_units": self.tryunits(ds, "coord3"),  # probably meters
                "coord1_units": self.tryunits(ds, "coord1"),  # probably meters
                "coord2_units": self.tryunits(ds, "coord2"),  # probably meters
            }
        else:
            raise ValueError("unknow type of dataset to make superdroplets from")

    def tryopen_dataset(self, dataset: Path):
        if isinstance(dataset, str) or isinstance(dataset, Path):
            print("supers dataset: ", dataset)
            return xr.open_dataset(dataset, engine="zarr", consolidated=False)
        else:
            return dataset

    def tryvar(self, ds: xr.Dataset, raggedcount: xr.DataArray, var: str):
        """attempts to return variable in form of ragged array
        (ak.Array) with dims [time, raggedcount]
        for a variable "var" in xarray dataset 'ds'.
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

    def get_units(self, var: str):
        return self._units[var + "_units"]

    def get_variable(self, var: str, units: bool):
        if not units:
            return self._raw_data[var]
        else:
            return self._raw_data[var], self.get_units(var)

    def sdId(self, units=False):
        return self.get_variable("sdId", units)

    def sdgbxindex(self, units=False):
        return self.get_variable("sdgbxindex", units)

    def xi(self, units=False):
        return self.get_variable("xi", units)

    def radius(self, units=False):
        return self.get_variable("radius", units)

    def msol(self, units=False):
        return self.get_variable("msol", units)

    def coord3(self, units=False):
        return self.get_variable("coord3", units)

    def coord1(self, units=False):
        return self.get_variable("coord1", units)

    def coord2(self, units=False):
        return self.get_variable("coord2", units)

    def time(self, units=False):
        if "time" in self._variables:
            return self.get_variable("time", units)
        print("time is not attached to superdrops")
        return None

    def radius_units(self):
        return self.get_units("radius")

    def msol_units(self):
        return self.get_units("msol")

    def coord3_units(self):
        return self.get_units("coord3")

    def coord1_units(self):
        return self.get_units("coord1")

    def coord2_units(self):
        return self.get_units("coord2")

    def time_units(self):
        if "time" in self._variables:
            return self.get_units("time")
        print("time is not attached to superdrops")
        return None

    def __getitem__(self, key):
        if key == "sdId":
            return self.sdId()
        elif key == "sdgbxindex":
            return self.sdgbxindex()
        elif key == "xi":
            return self.xi()
        elif key == "radius":
            return self.radius()
        elif key == "msol":
            return self.msol()
        elif key == "coord3":
            return self.coord3()
        elif key == "coord1":
            return self.coord1()
        elif key == "coord2":
            return self.coord2()
        elif key == "time":
            return self.time()
        elif key == "radius_units":
            return self.radius_units()
        elif key == "msol_units":
            return self.msol_units()
        elif key == "coord3_units":
            return self.coord3_units()
        elif key == "coord1_units":
            return self.coord1_units()
        elif key == "coord2_units":
            return self.coord2_units()
        elif key == "time_units":
            return self.time_units()
        else:
            err = "no known return provided for " + key + " key"
            raise ValueError(err)

    def attach_time(
        self, time: np.ndarray, time_units: str, do_reshape=False, var4reshape="sdId"
    ):
        if "time" not in self._raw_data:
            if do_reshape:
                # e.g. use sdId array to make 1-D time array have matching (ragged) dimension
                time = ak.broadcast_arrays(self._raw_data[var4reshape], time)[1]
            self._variables = tuple(list(self._variables) + ["time"])
            self._raw_data["time"] = ak.Array(time)
            self._units["time_units"] = time_units
        else:
            print("time already attached to superdrops")

    def detach_time(self):
        if "time" not in self._raw_data:
            print("time is not attached to superdrops")
        else:
            variables = list(self._variables)
            variables.remove("time")
            self._variables = tuple(variables)
            del self._raw_data["time"]
            del self._units["time_units"]

    def sample(self, sample_var: str, sample_values="all", variables2sample="all"):
        if isinstance(sample_var, str):
            sample_var = self._raw_data[sample_var]

        if isinstance(sample_values, str) and sample_values == "all":
            sample_values = list(np.unique(ak.flatten(sample_var)))
        elif not isinstance(sample_values, list):
            sample_values = [sample_values]

        if isinstance(variables2sample, str) and variables2sample == "all":
            variables2sample = self._variables
        elif not isinstance(variables2sample, list):
            variables2sample = [variables2sample]

        raw_data = {var: [] for var in variables2sample}
        for value in sample_values:
            mask = ak.Array(sample_var == value)
            for var in variables2sample:
                var_for_sample_value = ak.flatten(
                    ak.drop_none(self._raw_data[var].mask[mask])
                )
                raw_data[var].append(var_for_sample_value)

        superdrops_sample = Superdrops(self)
        superdrops_sample._variables = raw_data.keys()
        superdrops_sample._raw_data = raw_data

        return superdrops_sample

    def random_sample(
        self,
        sample_var: str,
        ndrops2sample: int,
        sample_population="all",
        variables2sample="all",
    ):
        if isinstance(sample_var, str):
            sample_var = self._raw_data[sample_var]

        if sample_population == "all":
            sample_population = list(np.unique(ak.flatten(sample_var)))

        sample_values = random.sample(sample_population, ndrops2sample)
        return self.sample(
            sample_var, sample_values=sample_values, variables2sample=variables2sample
        )

    def indexes_of_time_slices(self, time, times2select):
        if isinstance(time, list):
            time_indexes = []
            for i in range(len(time)):
                time_indexes.append(self.indexes_of_time_slices(time[i], times2select))
            return time_indexes

        if len(time) == ak.count(time):  # time is 1-D array
            time_indexes = []
            for t in times2select:
                idx = ak.argmin(abs((time - t)))
                time_indexes.append(idx)
        else:  # time assumed to be ragged array like superdrops data
            time_indexes = []
            for t in times2select:
                idx = ak.nanargmin(abs((time - t)), axis=0)[0]
                time_indexes.append(idx)

        return time_indexes

    def time_slice(
        self,
        times2select,
        variables2slice="all",
        attach_time=False,
        time=None,
        time_units=None,
    ):
        if attach_time:
            self.attach_time(time, time_units, do_reshape=False)

        if isinstance(variables2slice, str) and variables2slice == "all":
            variables2slice = self._variables
        elif not isinstance(variables2slice, list):
            variables2slice = [variables2slice]

        time_indexes = self.indexes_of_time_slices(self.time(), times2select)

        raw_data = {var: [] for var in variables2slice}
        for var in variables2slice:
            if isinstance(
                self._raw_data[var], list
            ):  # superdrops sample is being sliced
                time_slice_var = []
                for i in range(len(self._raw_data[var])):
                    time_slice_var.append(self._raw_data[var][i][time_indexes[i]])
            else:
                time_slice_var = self._raw_data[var][time_indexes]
            raw_data[var] = time_slice_var

        superdrops_timeslice = Superdrops(self)
        superdrops_timeslice._variables = raw_data.keys()
        superdrops_timeslice._raw_data = raw_data

        return superdrops_timeslice
