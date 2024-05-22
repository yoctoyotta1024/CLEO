"""
Copyright (c) 2024 MPI-M, Clara Bayley


----- CLEO -----
File: supersdata.py
Project: sdmout_src
Created Date: Tuesday 24th October 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
Last Modified: Tuesday 7th May 2024
Modified By: CB
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
python class to handle superdroplet attributes data from SDM zarr store in ragged array format
"""
# %%
import numpy as np
import xarray as xr
import awkward as ak
import warnings
import os
from typing import Union
import pytest


def akward_array_to_lagrange_array(
    data: ak.Array,
    dim1: ak.Array,
    dim2: ak.Array,
    dim1_as_index=False,
    check_indices_uniqueness=False,
):
    """
    This function converts a variable of the superdroplet dataset to a numpy array with the dimensions of the superdroplet dataset.
    The function assumes that the variable is a scalar value for each superdroplet at each time step.

    If you want to use it with the SupersData class, you can use the following syntax:
    >>> akward_array_to_lagrange_array(sddata[varname], sddata.time, sddata["sdId"])

    It will create a regular numpy array with the dimensions of the superdroplet dataset.
    N : number of superdroplets
    T : number of time steps
    (T,N) : shape of the output array


    The maximum "N" of ``dim2`` values is used to create the number of columns in the output array!
    The length "T" of ``dim1`` is used to create the number of columns in the output array!

    The output array the has shape (T,N).

    The input arrays are assumed to have the dimensions:
     - ``data``: T * var * dtype
     - ``dim1``: T * int
     - ``dim2``: T * var * int   (var is the number of superdroplets at each time step, which must be smaller than "N"!)

    Values in ``dim2`` need to be the indices for the columns of the output array!
    Values in ``dim1`` are not used as indices by default!
    If ``dim1_as_index`` is set to ``True``, the values of ``dim1`` are used as indices for the rows of the output array!


    Parameters
    ----------
    data : ak.Array
        The variable of the superdroplet dataset. The variable must be a scalar value for each superdroplet at each time step.
    dim1 : ak.Array
        The first dimension of the output numpy array. This is usually the time dimension.
        It is not used as index by default.
        This can be changed by setting ``dim1_as_index`` to ``True``.
    dim2 : ak.Array
        The second dimension of the output numpy array. This is usually the superdroplet dimension.
    dim1_as_index : bool, optional
        If set to ``True``, the values of ``dim1`` are used as indices for the rows of the output array. The default is ``False``.
    check_indices_uniqueness : bool, optional
        If set to ``True``, the uniqueness of the indices is checked. The default is ``False``.
        If the indices are not unique, a ValueError is raised.
        This check is very slow for large arrays! Better make sure the indices are unique before calling this function.

    Returns
    -------
    np.ndarray
        A numpy array with the dimensions of the superdroplet dataset. The variable is stored in the numpy array.

    Raises
    ------
    ValueError
        If the number of superdroplets and the number of values in the variable do not match.
    ValueError
        If the number of superdroplets and the number of values in the variable do not match for all time steps.
    ValueError
        Only if ``check_indices_uniqueness`` is set to ``True``. If the indice tuples are not unique.

    Examples
    --------

    The typical use, where time and superdroplet indices are given as arrays and "time" is not used as index.

    >>> data = ak.Array([
            [10, 20, 30],
            [],
            [40, 50],
        ])
    >>> time = ak.Array([
            10.0,
            20.0,
            45.0,
        ])
    >>> id = ak.Array([
            [0, 1, 2],
            [],
            [2, 3],
        ])
    >>> akward_array_to_lagrange_array(data = data, dim1 = time, dim2 = id, dim1_as_index = False, check_indices_uniqueness = False)
    ... [10., 20., 30., nan],
    ... [nan, nan, nan, nan],
    ... [nan, nan, 40., 50.]])

    In the following example, the time_index is used. As seen, the combination (0,0) from "time_index"and "id" is given twice. The last value will overwrite the previous.
    The check for uniqueness would through a ValueError.

    >>> data = ak.Array([
            [10, 20, 30],
            [],
            [40, 50],
            [90],
        ])
    >>> time_index = ak.Array([
            0,
            1,
            2,
            0,
        ])
    >>> id = ak.Array([
            [0, 1, 2],
            [],
            [2, 3],
            [0],
        ])
    >>> akward_array_to_lagrange_array(data = data, dim1 = time_index, dim2 = id, dim1_as_index = True, check_indices_uniqueness = False)
    ... [[90. 20. 30. nan]
    ...  [nan nan nan nan]
    ...  [nan nan 40. 50.]]
    >>> akward_array_to_lagrange_array(data = data, dim1 = time_index, dim2 = id, dim1_as_index = True, check_indices_uniqueness = True)
    ... ValueError: The indice tuples are not unique.
    ... This would lead to overwriting values in the numpy array.
    ... The reason might be, that the time indices aren't unique along axis 0 already.

    """

    # create the output dimensions of the numpy array which are necessary to store the data.

    if dim1_as_index is False:
        T = int(ak.num(dim1, axis=0))
        time_index = np.arange(T)
    else:
        time_index = dim1
        T = int(ak.max(dim1) + 1)
    # The superdroplets are identified by their id.
    # Use the maximum value of the superdroplet index
    # The ids start with id "0", so the maximum id is the number of superdroplets - 1!
    N = int(ak.max(dim2) + 1)
    superdroplet_index = dim2

    if ak.count(superdroplet_index) != ak.count(data):
        raise ValueError(
            f"The number of superdroplets ({ak.count(superdroplet_index)}) and the number of values in the variable ({ak.count(data)}) do not match"
        )
    if not ak.all(ak.num(data) == ak.num(data)):
        raise ValueError(
            "The number of superdroplets and the number of values in the variable do not match for all time steps"
        )

    # check if the resulting array is sparse and inform the User
    filled_percentage = ak.count(superdroplet_index) / (N * T) * 100
    # print(f"{filled_percentage:.2f} % of the regular array is filled with values. Total number of values is {ak.count(superdroplet_index)} out of {N * T} possible values.")
    if filled_percentage < 50:
        warnings.warn(
            "The resulting array is sparse. This might lead to significant memory usage"
        )

    # create tuples of all datapoint in the variable's akward array
    # For this, a cartesian product of the time_index and the superdroplet_index is created
    # The cartesian product is a tuple of all possible combinations of the two arrays
    # It is important to do this along axis 0 (time dimension). Otherwise, only unique combinations are created
    # The resulting array is then flattened to have a list of tuples
    # The list of tuples is then unzipped, to seperate the time and superdroplet indeices into two arrays
    i, j = ak.unzip(ak.flatten(ak.cartesian((time_index, superdroplet_index), axis=1)))

    # Check if the indice tuples are unique.
    # For this, one can simply compute the index value they represented in a raveled array.
    # so i * N + j should be unique for all i, j
    # flattened_index = i * N + j
    if check_indices_uniqueness is True:
        if not ak.count(i * N + j) == ak.count(np.unique(i * N + j)):
            raise ValueError(
                "The indice tuples are not unique.\n"
                + "This would lead to overwriting values in the numpy array.\n"
                + "The reason might be, that the time indices aren't unique along axis 0 already.\n"
            )
    else:
        warnings.warn(
            "The uniqueness of the indices is not checked. This might lead to overwriting values in the numpy array."
        )

    result_numpy = np.empty((T, N)) * np.nan
    result_numpy[i, j] = ak.flatten(data)
    return result_numpy


class SuperdropProperties:
    """Contains attributes common to all superdroplets and functions
    for calculating derived ones"""

    def __init__(self, consts):
        """Common attributes shared by superdroplets"""

        # density of liquid in droplets (=density of water at 300K) [Kg/m^3]
        self.RHO_L = consts["RHO_L"]

        # droplet solute properties
        self.RHO_SOL = consts["RHO_SOL"]  # density of (dry) solute [Kg/m^3]
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
    def __init__(self, dataset: os.PathLike, consts: dict):
        SuperdropProperties.__init__(self, consts)

        self.ds = self.tryopen_dataset(dataset)
        self.rgdcount = self.ds["raggedcount"].values  # ragged count variable

        self.time = ak.Array(self.ds.time.data)
        self.sdId = self.tryvar(self.ds, self.rgdcount, "sdId")
        self.sdgbxindex = self.tryvar(self.ds, self.rgdcount, "sdgbxindex")
        self.xi = self.tryvar(self.ds, self.rgdcount, "xi")
        self.radius = self.tryvar(self.ds, self.rgdcount, "radius")
        self.msol = self.tryvar(self.ds, self.rgdcount, "msol")

        self.coord3 = self.tryvar(self.ds, self.rgdcount, "coord3")
        self.coord1 = self.tryvar(self.ds, self.rgdcount, "coord1")
        self.coord2 = self.tryvar(self.ds, self.rgdcount, "coord2")

        # probably microns ie. 'micro m'
        self.radius_units = self.tryunits(self.ds, "radius")
        self.msol_units = self.tryunits(self.ds, "msol")  # probably gramms
        self.coord3_units = self.tryunits(self.ds, "coord3")  # probably meters
        self.coord1_units = self.tryunits(self.ds, "coord1")  # probably meters
        self.coord2_units = self.tryunits(self.ds, "coord2")  # probably meters

    def tryopen_dataset(self, dataset):
        if isinstance(dataset, str):
            print("supers dataset: ", dataset)
            return xr.open_dataset(dataset, engine="zarr", consolidated=False)
        else:
            return dataset

    def tryvar(self, ds, raggedcount, var):
        """attempts to return variable in form of ragged array
        (ak.Array) with dims [time, raggedcount]
        for a variable "var" in xarray dataset 'ds'.
        If attempt fails, returns empty array instead"""
        try:
            return ak.unflatten(ds[var].values, raggedcount)
        except KeyError:
            return ak.Array([])

    def tryunits(self, ds, var):
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
        elif key == "sdgbxindex_units":
            return self.sdgbxindex_units
        elif key == "xi":
            return self.xi
        elif key == "xi_units":
            return self.xi_units
        elif key == "radius":
            return self.radius
        elif key == "radius_units":
            return self.radius_units
        elif key == "msol":
            return self.msol
        elif key == "msol_units":
            return self.msol_units
        elif key == "coord3":
            return self.coord3
        elif key == "coord3_units":
            return self.coord3_units
        elif key == "coord1":
            return self.coord1
        elif key == "coord1_units":
            return self.coord1_units
        elif key == "coord2":
            return self.coord2
        elif key == "coord2_units":
            return self.coord2_units
        else:
            err = "no known return provided for " + key + " key"
            raise ValueError(err)

    def variable_to_regular_array(
        self,
        varname: str,
        dtype: type = np.ndarray,
        metadata: dict = dict(),
        check_indices_uniqueness: bool = False,
    ):
        """
        This function converts an awkward array to a regular array (either numpy ndarray or xarray DataArray)
        based on the provided dtype. The conversion is done for a specific variable name in the dataset.

        Parameters:
        varname (str): The name of the variable in the dataset to be converted.
        dtype (type, optional): The type of the output array. It can be either numpy ndarray or xarray DataArray.
                                Default is numpy ndarray.

        Returns:
        numpy.ndarray or xarray.DataArray: The converted regular array.

        Raises:
        ValueError: If the dtype provided is not supported. Only numpy ndarray and xarray DataArray are supported.
        """
        if dtype == np.ndarray:
            return akward_array_to_lagrange_array(
                self[varname],
                self.time,
                self["sdId"],
                check_indices_uniqueness=check_indices_uniqueness,
            )
        elif dtype == xr.DataArray:
            regular_array = self.variable_to_regular_array(
                varname=varname,
                dtype=np.ndarray,
                check_indices_uniqueness=check_indices_uniqueness,
            )
            result = xr.DataArray(
                regular_array,
                dims=["time", "sdId"],
                coords={
                    "time": self.time.to_list(),
                    "sdId": np.arange(regular_array.shape[1]),
                },
                name=varname,
                attrs=metadata,
            )
            return result
        else:
            raise ValueError("dtype is not supported. Use np.ndarray, xr.DataArray")

    def to_Dataset(self, check_indices_uniqueness: bool = False):
        """
        This function converts all variables in the dataset to regular arrays (either numpy ndarray or xarray DataArray)
        and combines them to a xarray Dataset.

        Returns
        -------
        xarray.Dataset
            The dataset with all variables converted to regular arrays.
            The dimensions are "time" and "sdId" for all variables.
            The output variables are "sdId", "sdgbxindex", "xi", "radius", "coord1", "coord2", "coord3".
        """
        varnames = ["sdId", "sdgbxindex", "xi", "radius", "coord1", "coord2", "coord3"]
        result_dict = dict()
        for varname in varnames:
            try:
                metadata = dict()
                # try to add units:
                try:
                    metadata["units"] = self[varname + "_units"]
                except ValueError:
                    print("No units found for", varname)
                    pass
                # create the regular array
                result_dict[varname] = self.variable_to_regular_array(
                    varname=varname,
                    dtype=xr.DataArray,
                    metadata=metadata,  # these are the metadata
                    check_indices_uniqueness=check_indices_uniqueness,  # check for uniqueness of indices
                )
            except ValueError:
                print(f"Could not create regular array for {varname}")
        # combine all variables to a dataset
        result = xr.Dataset(result_dict)
        result["time"].attrs = self.ds["time"].attrs
        result.attrs.update(**self.ds.attrs)
        result.attrs[
            "description"
        ] = "Regular shaped arrays created from the awkward arrays of the superdroplet dataset."
        result.attrs["creator_2"] = "Nils Niebaum for the regular array conversion."

        return result


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


IntegerArray = Union[np.ndarray, ak.highlevel.Array]


def get_awkward_shape(a):
    """
    Get the shape of the awkward array a as a list.
    Variable axis lengths are replaced by ``np.nan``.

    Parameters
    ----------
    a : ak.Array
        The input array.

    Returns
    -------
    list
        The shape of the array as a list.
        ``var`` positions are replaced by ``np.nan``.
    """

    # check for number of dimensions
    ndim = a.ndim
    # create output list
    shape = []
    # For each dinemsion, get the number of elements.
    # If the number of elements changes over the axis, np.nan is used to indicate a variable axis length
    for dim in range(ndim):
        num = ak.num(a, axis=dim)
        # for the 0th axis, the number of elements is an integer
        if isinstance(num, np.ndarray):
            num = int(num)
            shape.append(num)
        else:
            maxi = int(ak.max(num))
            mini = int(ak.min(num))
            if maxi == mini:
                shape.append(maxi)
            else:
                shape.append(np.nan)
    return shape


@pytest.mark.parametrize(
    "a, should",
    [
        (ak.Array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10]), [10]),
        (ak.Array([[1, 2, 3], [4, 5, 6]]), [2, 3]),
        (ak.Array([[[1, 2], [3, 4]], [[5, 6], [7, 8]]]), [2, 2, 2]),
        (ak.Array([[1, 2, 3], [4, 5, 6, 7]]), [2, np.nan]),
    ],
)
def test_get_awkward_shape(a, should):
    print(get_awkward_shape(a))
    assert get_awkward_shape(a) == should


def ak_digitize_2D(
    x: ak.highlevel.Array, bins: np.ndarray, right: bool = False
) -> ak.highlevel.Array:
    """
    This function takes a 2D awkward array and bins it into the provided bins.
    The binning is done using ``numpy.digitize`` along the flattened array.

    Note
    ----
    The input array must be 2D.
    None and np.nan will be stored in the last bin.

    Parameters
    ----------
    x : ak.highlevel.Array
        The input array to be binned.
    bins : np.ndarray
        The bins to be used for binning the array.
    right : bool, optional
        As from the numpy documentation:
        Indicating whether the intervals include the right or the left bin edge.
        Default behavior is (right==False) indicating that the interval does not include the right edge.
        The left bin end is open in this case, i.e., bins[i-1] <= x < bins[i] is the default behavior for monotonically increasing bins.

    Returns
    -------
    indices : ak.highlevel.Array of ints
        Output array of indices, of same shape as x.

    Example
    -------
    >>> x = ak.Array([[1, 2, 3], [4, 5, 6]])
    >>> bins = np.array([0, 2, 4, 6])
    >>> array_to_bin_index(x, bins, right = False)
    [[1, 2, 2],
     [3, 3, 4]]
    ---------------------
    type: 2 * var * int64

    References
    ----------
    https://numpy.org/doc/stable/reference/generated/numpy.digitize.html
    """

    if x.ndim != 2:
        raise ValueError("Array x must be 2D")
    digi, num = ak.flatten(x), ak.num(x)

    if bins.ndim != 1:
        raise ValueError("Bins must be 1D")

    digi = np.digitize(x=digi, bins=bins, right=right)
    digi = ak.unflatten(digi, num)
    return digi


@pytest.mark.parametrize(
    "x, bins, should",
    [
        (
            ak.Array([[1, 2, 3, -1, 101], [4, 5, 6, np.nan, 90, None]]),
            np.array([0, 3, 6, 100]),
            ak.Array([[1, 1, 2, 0, 4], [2, 2, 3, 4, 3, 4]]),
        ),
    ],
)
def test_ak_digitize_2D(x, bins, should):
    """Tests returns of the ak_digitize_2D function"""
    x_binned = ak_digitize_2D(x, bins)
    # check for equality
    assert ak.sum(x_binned != should) == 0
    # check for same shape
    assert get_awkward_shape(x_binned) == get_awkward_shape(should)
    # check for same counts
    assert ak.count(x_binned) == ak.count(should)


@pytest.mark.parametrize(
    "x, bins, exception",
    [
        (
            ak.Array([[[1, 2, 3, -1, 101], [4, 5, 6, np.nan, 90, None]], [None]]),
            np.array([0, 3, 6, 100]),
            ValueError,
        ),
        (
            ak.Array([[1, 2, 3, -1, 101], [4, 5, 6, np.nan, 90, None]]),
            np.array([[0], [3]]),
            ValueError,
        ),
    ],
)
def test_array_to_bin_index_exception(x, bins, exception):
    """Tests for exceptions in the ak_digitize_2D function"""
    with pytest.raises(exception):
        ak_digitize_2D(x, bins)


def ak_digitize_3D(
    x: ak.highlevel.Array, bins: np.ndarray, right: bool = False
) -> ak.highlevel.Array:
    """
    This function takes a 3D awkward array and bins it into the provided bins.
    The binning is done using ``numpy.digitize`` along the flattened array.

    Note
    ----
    The input array must be 3D.
    None and np.nan will be stored in the last bin.
    Only values in the last axis are allowed! This ``ak.Array([[[1, 2, 3, -1, 101], [4, 5, 6, np.nan, 90, None]], [None, 1]])`` is not valid!


    Parameters
    ----------
    x : ak.highlevel.Array
        The input array to be binned.
    bins : np.ndarray
        The bins to be used for binning the array.
    right : bool, optional
        As from the numpy documentation:
        Indicating whether the intervals include the right or the left bin edge.
        Default behavior is (right==False) indicating that the interval does not include the right edge.
        The left bin end is open in this case, i.e., bins[i-1] <= x < bins[i] is the default behavior for monotonically increasing bins.

    Returns
    -------
    indices : ak.highlevel.Array of ints
        Output array of indices, of same shape as x.

    Example
    -------
    >>> x = ak.Array([[[1, 2, 3], [4, 5, 6]], [[1], [2]]])
    >>> bins = np.array([0, 2, 4, 6])
    >>> array_to_bin_index(x, bins, right = False)
    [
        [[1, 2, 2],[3, 3, 4]]
        [[1], [2]],
    ]
    ---------------------
    type: 2 * 2 * var * int64

    References
    ----------
    https://numpy.org/doc/stable/reference/generated/numpy.digitize.html
    """

    if x.ndim != 3:
        raise ValueError("Array x must be 3D")
    digi, num0 = ak.flatten(x), ak.num(x)
    digi, num1 = ak.flatten(digi), ak.num(digi)

    if bins.ndim != 1:
        raise ValueError("Bins must be 1D")

    digi = np.digitize(x=digi, bins=bins, right=right)
    digi = ak.unflatten(digi, num1)
    digi = ak.unflatten(digi, num0)
    return digi


@pytest.mark.parametrize(
    "x, bins, should",
    [
        (
            ak.Array([[[1, 2, 3], [4, 5]], [[1], [2]]]),
            np.array([0, 2, 4, 6]),
            ak.Array([[[1, 2, 2], [3, 3]], [[1], [2]]]),
        ),
    ],
)
def test_ak_digitize_3D(x, bins, should):
    """Tests returns of the ak_digitize_3D function"""
    x_binned = ak_digitize_3D(x, bins)
    # check for equality
    assert ak.sum(x_binned != should) == 0
    # check for same shape
    assert get_awkward_shape(x_binned) == get_awkward_shape(should)
    # check for same counts
    assert ak.count(x_binned) == ak.count(should)


@pytest.mark.parametrize(
    "x, bins, exception",
    [
        # The first test is not a valid case, because only values in the last axis are allowed
        (
            ak.Array([[[1, 2, 3, -1, 101], [4, 5, 6, np.nan, 90, None]], [None, 1]]),
            np.array([0, 3, 6, 100]),
            ValueError,
        ),
        (
            ak.Array([[1, 2, 3, -1, 101], [4, 5, 6, np.nan, 90, None, 1]]),
            np.array([[0], [3]]),
            ValueError,
        ),
    ],
)
def test_ak_digitize_3D_exception(x, bins, exception):
    """Tests for exceptions in the ak_digitize_3D function"""
    with pytest.raises(exception):
        ak_digitize_3D(x, bins)


class SupersAttribute:
    """
    This class is used to store the attributes of the superdroplets.
    It can be an indexer or a variable.

    The class has the following attributes:
    - name (str): The name of the attribute.
    - units (str): The units of the attribute.
    - data (ak.Array): The data of the attribute.
    - metadata (dict): The metadata of the attribute.
    """

    def __init__(
        self,
        name: str,
        data: Union[xr.Dataset, xr.DataArray, np.ndarray, ak.Array],
        units: str = "",
        metadata: dict = dict(),
    ):
        """
        The constructor of the class uses the data
        This function extracts an attribute from a dataset and stores it in a SupersAttribute object.

        Parameters
        ----------
        name : str
            The name of the attribute.
        data : xr.Dataset or xr.DataArray, or np.ndarray or ak.Array
            The data of the attribute.
            If a xr.Dataset or xr.DataArray is provided, the data is extracted from the dataset using the name of the attribute.
        units : str, optional
            The units of the attribute.
            Default is an empty string.
        metadata : dict, optional
            The metadata of the attribute.
            Default is an empty dictionary.
        """
        # store the name
        self.set_name(name)

        if isinstance(data, xr.Dataset):
            self.set_from_Dataset(data)
        elif isinstance(data, xr.DataArray):
            self.set_from_DataArray(data)
        else:
            self.set_data(data)
            self.set_units(units)
            self.set_metadata(metadata)

    def set_name(self, name: str):
        """
        This function sets the name of the attribute.
        The name is stored in the attribute name.

        Parameters
        ----------
        name : str
            The name of the attribute.
        """
        self.name = name

    def set_units(self, units):
        """
        This function sets the units of the attribute.
        The units are stored in the attribute units.

        Parameters
        ----------
        units : str
            The units of the attribute.
            If units are not provided, it is tried to fill the units from the dataset.
        """
        self.units = units

    def set_data(self, data: Union[ak.Array, np.ndarray]):
        """
        This function sets the data of the attribute.
        The data is stored in the attribute data.

        Parameters
        ----------
        data : ak.Array or np.ndarray or xr.Dataset
            The data of the attribute.
            If a xr.Dataset is provided, the data is extracted from the dataset using the name of the attribute.
        """
        if isinstance(data, np.ndarray):
            self.data = ak.Array(data)
        elif isinstance(data, ak.Array):
            self.data = data
        else:
            raise ValueError("Data must be an ak.Array or np.ndarray")

    def set_metadata(self, metadata: dict):
        """
        This function sets the metadata of the attribute.
        The metadata is stored in the attribute metadata.

        Parameters
        ----------
        metadata : dict
            The metadata of the attribute.
        """
        self.metadata = metadata

    def set_from_Dataset(self, ds: xr.Dataset):
        """
        This function sets the data of the attribute from a dataset.
        The data is stored in the attribute data.

        Parameters
        ----------
        ds : xr.Dataset
            The dataset containing the attribute.
        """
        self.set_data(data=ds[self.name].values)
        try:
            self.set_units(units=ds[self.name].units)
        except AttributeError:
            self.set_units(units="")
        try:
            self.set_metadata(metadata=ds[self.name].attrs)
        except AttributeError:
            self.set_metadata(metadata=dict())

    def set_from_DataArray(self, da: xr.DataArray):
        """
        This function sets the data of the attribute from a DataArray.
        The data is stored in the attribute data.

        Parameters
        ----------
        da : xr.DataArray
            The DataArray containing the attribute.
        """
        self.set_data(data=da.values)
        try:
            self.set_units(units=da.units)
        except AttributeError:
            self.set_units(units="")
        try:
            self.set_metadata(metadata=da.attrs)
        except AttributeError:
            self.set_metadata(metadata=dict())

    def __str__(self):
        """
        This function returns a string representation of the attribute.
        The string representation contains the name of the attribute and the units.

        Returns
        -------
        str
            The string representation of the attribute.
        """
        return f"{self.name} ({self.units})\n{self.data}"

    @classmethod
    def attribute_to_indexer(cls, attribute: "SupersAttribute") -> "SupersIndexer":
        """
        This function converts an attribute to an indexer.
        The attribute is converted to an indexer by creating a SupersIndexer object.

        Parameters
        ----------
        attribute : SupersAttribute
            The attribute to be converted to an indexer.

        Returns
        -------
        SupersIndexer
            The indexer created from the attribute.
        """
        return SupersIndexer(
            name=attribute.name,
            data=attribute.data,
            units=attribute.units,
            metadata=attribute.metadata,
        )

    @classmethod
    def attribute_to_indexer_binned(
        cls, attribute: "SupersAttribute", bin_edges: np.ndarray, right: bool = False
    ) -> "SupersIndexerBinned":
        """
        This function converts an attribute to a binned indexer.
        The attribute is converted to a binned indexer by creating a SupersIndexerBinned object.

        Parameters
        ----------
        attribute : SupersAttribute
            The attribute to be converted to a binned indexer.
        bin_edges : np.ndarray
            The bin edges of the binned indexer.
        right : bool, optional
            As from the numpy documentation:
            Indicating whether the intervals include the right or the left bin edge.
            Default behavior is (right==False) indicating that the interval does not include the right edge.

        Returns
        -------
        SupersIndexerBinned
            The binned indexer created from the attribute.
        """
        return SupersIndexerBinned(
            name=attribute.name,
            data=attribute.data,
            bin_edges=bin_edges,
            right=right,
            units=attribute.units,
            metadata=attribute.metadata,
        )


class SupersIndexer(SupersAttribute):
    """
    This class is used to store an indexers attribute of the superdroplets dataset.
    It is a subclass of the SupersAttribute class.
    This class can be used to store the indexer of the superdroplets dataset.

    It has the following attributes:
    - name (str): The name of the indexer.
    - units (str): The units of the indexer.
    - data (ak.Array): The data of the indexer.
    - digitized_data (ak.Array): The data of the indexer as digitized values.
    - metadata (dict): The metadata of the indexer.
    """

    def __init__(
        self,
        name: str,
        data: Union[xr.Dataset, np.ndarray, ak.Array],
        units: str = "",
        metadata: dict = dict(),
    ):
        """
        The constructor of the class uses the data
        This function extracts an attribute from a dataset and stores it in a SupersAttribute object.

        Parameters
        ----------
        name : str
            The name of the attribute.
        data : xr.Dataset or np.ndarray or ak.Array
            The data of the attribute.
            If a xr.Dataset is provided, the data is extracted from the dataset using the name of the attribute.
        units : str, optional
            The units of the attribute.
            Default is an empty string.
        metadata : dict, optional
            The metadata of the attribute.
            Default is an empty dictionary.
        """

        super().__init__(name=name, data=data, units=units, metadata=metadata)

        self.set_digitized_data()

    def set_digitized_data(self):
        """
        This function sets the digitized data of the indexer.
        The digitized data is stored in the attribute digitized_data.
        In this class, the digitzed data is the same as the data.
        So the indexer should be integer values.
        """

        # digitize the data
        self.digitized_data = self.data

    @classmethod
    def indexer_to_indexer_binned(
        cls, attribute: "SupersAttribute", bin_edges: np.ndarray, right: bool = False
    ) -> "SupersIndexerBinned":
        """
        This function converts an indexer to a binned indexer.
        The indexer is converted to a binned indexer by creating a SupersIndexerBinned object.

        Parameters
        ----------
        attribute : SupersAttribute
            The attribute to be converted to a binned indexer.
        bin_edges : np.ndarray
            The bin edges of the binned indexer.
        right : bool, optional
            As from the numpy documentation:
            Indicating whether the intervals include the right or the left bin edge.
            Default behavior is (right==False) indicating that the interval does not include the right edge.

        Returns
        -------
        SupersIndexerBinned
            The binned indexer created from the indexer.
        """
        return SupersIndexerBinned(
            name=attribute.name,
            data=attribute.data,
            bin_edges=bin_edges,
            right=right,
            units=attribute.units,
            metadata=attribute.metadata,
        )

    @classmethod
    def indexer_to_attribute(cls, indexer: "SupersIndexer") -> "SupersAttribute":
        """
        This function converts an indexer to an attribute.
        The indexer is converted to an attribute by creating a SupersAttribute object.

        Parameters
        ----------
        indexer : SupersIndexer
            The indexer to be converted to an attribute.

        Returns
        -------
        SupersAttribute
            The attribute created from the indexer.
        """
        return SupersAttribute(
            name=indexer.name,
            data=indexer.data,
            units=indexer.units,
            metadata=indexer.metadata,
        )


class SupersIndexerBinned(SupersIndexer):
    """
    This class is used to store a binned indexers attribute of the superdroplets dataset.
    It is a subclass of the SupersAttribute class.
    This class can be used to store the binned indexer of the superdroplets dataset.

    It has the following attributes:
    - name (str): The name of the binned indexer.
    - units (str): The units of the binned indexer.
    - data (ak.Array): The data of the binned indexer.
    - digitized_data (ak.Array): The data of the binned indexer as digitized values.
    - metadata (dict): The metadata of the binned indexer.
    - bin_edges (np.ndarray): The bin edges of the binned indexer.
    - bin_centers (np.ndarray): The bin centers of the binned indexer.
    """

    def __init__(
        self,
        name: str,
        data: Union[xr.Dataset, np.ndarray, ak.Array],
        bin_edges: np.ndarray,
        right: bool = False,
        units: str = "",
        metadata: dict = dict(),
    ):
        """
        The constructor of the class uses the data
        This function extracts an attribute from a dataset and stores it in a SupersAttribute object.

        Parameters
        ----------
        name : str
            The name of the attribute.
        data : xr.Dataset or np.ndarray or ak.Array
            The data of the attribute.
            If a xr.Dataset is provided, the data is extracted from the dataset using the name of the attribute.
        bin_edges : np.ndarray
            The bin edges of the binned indexer.
        right : bool, optional
            As from the numpy documentation:
            Indicating whether the intervals include the right or the left bin edge.
            Default behavior is (right==False) indicating that the interval does not include the right edge.
        units : str, optional
            The units of the attribute.
            Default is an empty string.
        metadata : dict, optional
            The metadata of the attribute.
            Default is an empty dictionary.
        """

        super().__init__(name=name, data=data, units=units, metadata=metadata)

        self.set_bins(bin_edges=bin_edges)
        self.set_digitized_data(bins=self.bin_edges, right=right)

    def set_bins(self, bin_edges: np.ndarray):
        """
        This function sets the bins of the binned indexer.
        The bins are stored in the attribute bins.
        The bin centers are stored in the attribute bin_centers.

        Parameters
        ----------
        bin_edges : np.ndarray
            The bin edges of the binned indexer.
        """

        # set the bins
        self.bin_edges = bin_edges
        # set the bin centers
        self.bin_centers = (self.bin_edges[:-1] + self.bin_edges[1:]) / 2

    def set_digitized_data(self, bins: np.ndarray, right: bool = False):
        """
        Sets a digitized version of the data.
        It uses the numpy digitize function to digitize the data.

        The digitized data is stored in the attribute digitized_data.

        Note

        Parameters
        ----------
        bins : np.ndarray
            The bins to digitize the data to.
        right : bool, optional
            As from the numpy documentation:
            Indicating whether the intervals include the right or the left bin edge.
            Default behavior is (right==False) indicating that the interval does not include the right edge.

        """

        # digitize the data
        if self.data.ndim == 1:
            self.digitized_data = np.digitize(x=self.data, bins=bins, right=right)
        elif self.data.ndim == 2:
            self.digitized_data = ak_digitize_2D(x=self.data, bins=bins, right=right)
        elif self.data.ndim == 3:
            self.digitized_data = ak_digitize_3D(x=self.data, bins=bins, right=right)
        else:
            raise NotImplementedError(
                "Only 1D, 2D and 3D arrays are supported for digitization till now."
            )


# %%


class SupersDataIndexed(SuperdropProperties):
    attribute_names = [
        "sdId",
        "sdgbxindex",
        "xi",
        "radius",
        "msol",
        "coord3",
        "coord1",
        "coord2",
    ]

    def __init__(self, dataset: Union[os.PathLike, xr.Dataset], consts: dict):
        """
        The constructor of the class uses the data
        This function extracts the attributes of the superdroplets dataset and stores them in a SupersDataIndexed object.

        """

        SuperdropProperties.__init__(self, consts=consts)

        self.ds = self.tryopen_dataset(dataset)
        self.raggedcount = self.ds["raggedcount"].values  # ragged count variable

        self.attributes = dict()
        self.set_attributes(attribute_names=self.attribute_names)

    def tryopen_dataset(self, dataset: Union[os.PathLike, xr.Dataset]) -> xr.Dataset:
        if isinstance(dataset, str):
            print("supers dataset: ", dataset)
            return xr.open_dataset(dataset, engine="zarr", consolidated=False)
        elif isinstance(dataset, xr.Dataset):
            return dataset
        else:
            raise ValueError("dataset must be a path or a xarray dataset")

    def set_attributes(self, attribute_names):
        """
        This function sets the attributes of the superdroplets dataset.
        The attributes are stored in the attribute attributes.

        Parameters
        ----------
        attribute_names : list
            The names of the attributes to be stored.
        """

        # set all attribtes
        for name in attribute_names:
            self.set_attribute(name=name)

        # set time attribute
        self.set_time()

    def set_attribute(self, name: str):
        """
        This function sets an attribute of the superdroplets dataset.
        The attribute is stored in the attribute attributes.

        Parameters
        ----------
        name : str
            The name of the attribute to be stored.
        """

        self.attributes[name] = SupersAttribute(name=name, data=self.ds[name])

    def set_time(self):
        """
        This function sets the time attribute of the superdroplets dataset.
        The time attribute is stored in the attribute time.
        """

        self.time = np.repeat(
            ak.Array(self.ds.time.data),
            self.raggedcount,
        )

    def __getitem__(self, key):
        try:
            self.attributes[key]
        except KeyError:
            err = "no known return provided for " + key + " key"
            raise ValueError(err)

    def __str__(self):
        """
        This function returns a string representation of the superdroplets dataset.
        The string representation contains the names of the attributes and their units.

        Returns
        -------
        str
            The string representation of the superdroplets dataset.
        """
        result = ""
        for name, attribute in self.attributes.items():
            result += str(attribute) + "\n"
        return result
