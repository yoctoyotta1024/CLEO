'''
----- CLEO -----
File: supersdata.py
Project: sdmout_src
Created Date: Tuesday 24th October 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
Last Modified: Monday 20th November 2023
Modified By: CB
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
Copyright (c) 2023 MPI-M, Clara Bayley
-----
File Description:
python class to handle superdroplet
attributes data from SDM zarr store
in ragged array format
'''

import numpy as np
import xarray as xr
import awkward as ak
import warnings

def akward_array_to_lagrange_array(data : ak.Array, dim1 : ak.Array, dim2 : ak.Array, dim1_as_index = False, check_indices_uniqueness = False) :
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
        T = int(ak.num(dim1, axis = 0))
        time_index = np.arange(T)
    else :
        time_index = dim1
        T = int(ak.max(dim1) + 1)
    # The superdroplets are identified by their id.
    # Use the maximum value of the superdroplet index
    # The ids start with id "0", so the maximum id is the number of superdroplets - 1!
    N = int(ak.max(dim2) + 1)
    superdroplet_index = dim2

    if ak.count(superdroplet_index) != ak.count(data):
        raise ValueError(f"The number of superdroplets ({ak.count(superdroplet_index)}) and the number of values in the variable ({ak.count(data)}) do not match")
    if not ak.all(ak.num(data) == ak.num(data)):
        raise ValueError(f"The number of superdroplets and the number of values in the variable do not match for all time steps")

    # check if the resulting array is sparse and inform the User
    filled_percentage = ak.count(superdroplet_index) / (N * T) * 100
    # print(f"{filled_percentage:.2f} % of the regular array is filled with values. Total number of values is {ak.count(superdroplet_index)} out of {N * T} possible values.")
    if filled_percentage < 50:
        warnings.warn(f"The resulting array is sparse. This might lead to significant memory usage")

    # create tuples of all datapoint in the variable's akward array
    # For this, a cartesian product of the time_index and the superdroplet_index is created
    # The cartesian product is a tuple of all possible combinations of the two arrays
    # It is important to do this along axis 0 (time dimension). Otherwise, only unique combinations are created
    # The resulting array is then flattened to have a list of tuples
    # The list of tuples is then unzipped, to seperate the time and superdroplet indeices into two arrays
    i, j = ak.unzip(ak.flatten(ak.cartesian((time_index, superdroplet_index), axis = 1)))


    # Check if the indice tuples are unique.
    # For this, one can simply compute the index value they represented in a raveled array.
    # so i * N + j should be unique for all i, j
    # flattened_index = i * N + j
    if check_indices_uniqueness is True :
        if not ak.count(i * N + j) == ak.count(np.unique(i * N + j)) :
            raise ValueError(
                "The indice tuples are not unique.\n"\
                +"This would lead to overwriting values in the numpy array.\n"\
                +"The reason might be, that the time indices aren't unique along axis 0 already.\n"
                )
    else :
        warnings.warn("The uniqueness of the indices is not checked. This might lead to overwriting values in the numpy array.")

    result_numpy = np.empty((T, N)) * np.nan
    result_numpy[i, j] = ak.flatten(data)
    return result_numpy



class SuperdropProperties():
    '''Contains attributes common to all superdroplets and functions
    for calculating derived ones'''

    def __init__(self, consts):
        '''Common attributes shared by superdroplets'''

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
        ''' calculates effective density [g m^-3] of
      droplet such that mass_droplet, m = 4/3*pi*r^3 * rhoeff
      taking into account mass of liquid and mass of
      solute assuming solute occupies volume it
      would given its (dry) density, RHO_SOL. '''

        msol = msol/1000 # convert from grams to Kg
        r = r/1e6 # convert microns to m

        solfactor = 3*msol/(4.0*np.pi*(r**3))
        rhoeff = self.RHO_L + solfactor*(1-self.RHO_L/self.RHO_SOL)

        return rhoeff * 1000 #[g/m^3]

    def vol(self, r):
        ''' volume of droplet [m^3] '''

        r = r/1e6 # convert microns to m

        return 4.0/3.0 * np.pi * r**3

    def mass(self, r, msol):
        '''
        total mass of droplet (water + (dry) areosol) [g],
        m =  4/3*pi*rho_l**3 + msol(1-rho_l/rho_sol)
        ie. m = 4/3*pi*rhoeff*R**3
        '''

        msol = msol/1000 # convert from grams to Kg
        r = r/1e6 # convert microns to m

        msoleff = msol*(1-self.RHO_L/self.RHO_SOL) # effect of solute on mass
        m = msoleff + 4/3.0*np.pi*(r**3)*self.RHO_L

        return m * 1000 # [g]

    def m_water(self, r, msol):
        ''' mass of only water in droplet [g]'''

        msol = msol/1000 # convert msol from grams to Kg
        r = r/1e6 # convert microns to m

        v_sol = msol/self.RHO_SOL
        v_w = 4/3.0*np.pi*(r**3) - v_sol

        return self.RHO_L*v_w * 1000 #[g]

class SupersData(SuperdropProperties):

    def __init__(self, dataset, consts):
        SuperdropProperties.__init__(self, consts)

        self.ds = self.tryopen_dataset(dataset)
        self.rgdcount = self.ds["rgd_totnsupers"].values  # ragged count variable

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

        if type(dataset) == str:
            print("supers dataset: ", dataset)
            return xr.open_dataset(dataset, engine="zarr", consolidated=False)
        else:
            return dataset

    def tryvar(self, ds, raggedcount, var):
        ''' attempts to return variable in form of ragged array
        (ak.Array) with dims [dim_temp, raggedcount]
        for a variable "var" in xarray dataset 'ds'.
        If attempt fails, returns empty array instead '''
        try:
            return ak.unflatten(ds[var].values, raggedcount)
        except:
            return ak.Array([])

    def tryunits(self, ds, var):
        ''' attempts to return the units of a variable
        in xarray dataset 'ds'. If attempt fails, returns null '''
        try:
            return ds[var].units
        except:
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
            err = "no known return provided for "+key+" key"
            raise ValueError(err)

    def variable_to_regular_array(self, varname : str, dtype = np.ndarray) :
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
            return akward_array_to_lagrange_array(self[varname], self.time, self["sdId"])
        elif dtype == xr.DataArray:
            regular_array = self.variable_to_regular_array(varname=varname, dtype=np.ndarray)
            result = xr.DataArray(
                regular_array,
                dims = ["time", "sdId"],
                coords = {"time" : self.time.to_list(), "sdId" : np.arange(regular_array.shape[1])},
                name = varname
            )
            return result
        else:
            raise ValueError("dtype is not supported. Use np.ndarray, xr.DataArray")

    def to_Dataset(self) :
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
        result_list = []
        for varname in varnames:
            try :
                result_list.append(self.variable_to_regular_array(
                    varname = varname,
                    dtype = xr.DataArray
                    )
                )
            except ValueError:
                print(f"Could not create regular array for {varname}")
        result = xr.Dataset(dict(zip(varnames, result_list)))
        return result
class RainSupers(SuperdropProperties):

    def __init__(self, sddata, rlim=40):
        ''' return data for (rain)drops with radii > rlim.
        Default minimum raindrops size is rlim=40microns'''

        if type(sddata) != SupersData:
            sddata = SupersData(dataset=sddata)

        israin = sddata.radius >= rlim  # ak array True for raindrops

        self.totnsupers_rain = ak.num(israin[israin == True])
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
            err = "no known return provided for "+key+" key"
            raise ValueError(err)
