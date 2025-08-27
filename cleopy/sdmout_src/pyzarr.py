"""
Copyright (c) 2025 MPI-M, Clara Bayley


----- CLEO -----
File: pyzarr.py
Project: sdmout_src
Created Date: Tuesday 24th October 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
functions to return zarr data in useful formats for plotting e.g. for handling ragged arrays
of superdroplet attributes
"""


import numpy as np
import xarray as xr
import awkward as ak
from pathlib import Path

from . import thermodata
from . import superdrops
from . import massmoms
from . import timedata


def get_rawdataset(dataset):
    return xr.open_dataset(dataset, engine="zarr", consolidated=False)


def get_rawdata4key(dataset, key):
    return get_rawdataset(dataset)[key]


def get_rawdata4raggedkey(dataset, key):
    ds = get_rawdataset(dataset)
    return ak.unflatten(ds[key].values, ds["raggedcount"].values)


def raggedvar_fromzarr(ds, raggedcount, var):
    """returns ragged ak.Array dims [time, ragged]
    for a variable "var" in zarr ds"""
    if isinstance(ds, str) or isinstance(ds, Path):
        ds = get_rawdataset(ds)

    return ak.unflatten(ds[var].values, raggedcount)


def var4d_fromzarr(ds, ntime, ndims, key):
    """' returns 4D variable with dims
    [time, y, x, z] from zarr dataset "ds" """
    if isinstance(ds, str) or isinstance(ds, Path):
        ds = get_rawdataset(ds)

    reshape = [ntime] + list(ndims)
    return np.reshape(ds[key].values, reshape)


def var3d_fromzarr(ds, ndims, key):
    """' returns 3D variable with dims
    [y, x, z] from zarr dataset "ds" """
    if isinstance(ds, str) or isinstance(ds, Path):
        ds = get_rawdataset(ds)

    return np.reshape(ds[key].values, ndims)


def get_thermodata(dataset, ntime, ndims, consts, getwinds=False):
    """returns a thermodynamic data in a dictionary. The value under
    each key is the thermodynamics data in a 2D array
    with dimensions [time, gridbox]. E.g. thermo["qvap"][:,0] gives the
    timeseries of qvap for the 0th gridbox. thermo["qvap][0] gives
    the qvap of all gridboxes at the 0th output time"""

    thermo = thermodata.Thermodata(dataset, ntime, ndims, consts)
    if getwinds:
        winds = thermodata.Winddata(dataset, ntime, ndims, consts)
        return thermo, winds
    else:
        return thermo


def get_supers(dataset, consts, is_print=False):
    return superdrops.Superdrops(dataset, consts, is_print=is_print)


def get_time(dataset):
    return timedata.Time(dataset)


def get_massmoms(dataset, ntime, ndims):
    return massmoms.MassMoms(dataset, ntime, ndims)


def get_rainmassmoms(dataset, ntime, ndims):
    return massmoms.MassMoms(dataset, ntime, ndims, lab="_raindrops")


def get_gbxindex(dataset, ndims):
    return var3d_fromzarr(dataset, ndims, "gbxindex")


def get_totnsupers(dataset):
    if isinstance(dataset, str) or isinstance(dataset, Path):
        dataset = get_rawdataset(dataset)
    try:
        return dataset["totnsupers"].values
    except KeyError:
        return dataset["raggedcount"].values


def get_nsupers(dataset, ntime, ndims):
    return var4d_fromzarr(dataset, ntime, ndims, "nsupers")


def surfaceprecip_estimate(dataset, gbxs):
    """use last radius of SDs before they leave the domain to
    estimate the volume of precipitation at each timestep.
    Values should be approx. equal to sum over gbxs (multiplied
    by area_gbx/area_domain) of logbook values for precip"""

    ds = get_rawdataset(dataset)

    sdId = ak.unflatten(ds["sdId"].values, ds["raggedcount"].values)
    radius = ak.unflatten(ds["radius"].values, ds["raggedcount"].values)
    xi = ak.unflatten(ds["xi"].values, ds["raggedcount"].values)

    r3sum = []
    for ti in range(ds.time.shape[0] - 1):
        sd_ti, r_ti, xi_ti = sdId[ti], radius[ti], xi[ti]
        sds_gone = set(sd_ti) - set(
            sdId[ti + 1]
        )  # set of SD indexes that have left domain during timestep ti -> ti+1
        isgone = np.where(
            np.isin(sd_ti, list(sds_gone))
        )  # indexes in ragged arrays of SDs that leave during timestep ti -> ti+1
        r3sum.append(
            np.dot(r_ti[isgone] ** 3, xi_ti[isgone])
        )  # sum of (real) droplet radii^3 that left domain [microns^3]
    precipvol = (
        4 / 3 * np.pi * np.asarray(r3sum) / (1e18)
    )  # volume of water that left domain [m^3]

    domainy = np.amax(gbxs["yhalf"]) - np.amin(gbxs["yhalf"])  # [m]
    domainx = np.amax(gbxs["xhalf"]) - np.amin(gbxs["xhalf"])  # [m]
    deltat = np.diff(ds["time"].values) / 60 / 60  # [hrs]
    preciprate = precipvol * 1000 / (domainx * domainy) / deltat  # [mm/hr]

    precipaccum = np.cumsum(preciprate * deltat)  # [mm]
    preciprate = np.insert(preciprate, 0, 0)  # at t=0, precip rate = 0
    precipaccum = np.insert(precipaccum, 0, 0)  # at t=0, accumulated precip = 0

    return preciprate, precipaccum  # [mm/hr] , [mm]
