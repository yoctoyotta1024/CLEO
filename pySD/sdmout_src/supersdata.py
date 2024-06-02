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
# import libraries

import numpy as np
import xarray as xr
import awkward as ak
import os
from typing import Union, Tuple, Callable

from pySD.sdmout_src import sdtracing


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
            return sdtracing.akward_array_to_lagrange_array(
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

    def get_data(self):
        """
        This function returns the data of the attribute.

        Returns
        -------
        ak.Array
            The data of the attribute.
        """
        return self.data

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
        return f"{self.name} ({self.units})\n{self.data.typestr}"

    def __add__(self, other):
        """
        This method overloads the + operator.
        It performs element-wise addition of the data attribute with another object.

        Parameters
        ----------
        other : object
            The object to be added to the data attribute.

        Returns
        -------
        SupersAttribute
            A new SupersAttribute object with the result of the addition.
        """
        if isinstance(other, SupersAttribute):
            new_data = self.data + other.data
        else:
            new_data = self.data + other
        return SupersAttribute(
            name=self.name,
            data=new_data,
            units=self.units,
            metadata=self.metadata,
        )

    def __sub__(self, other):
        """
        This method overloads the - operator.
        It performs element-wise subtraction of the data attribute with another object.

        Parameters
        ----------
        other : object
            The object to be subtracted from the data attribute.

        Returns
        -------
        SupersAttribute
            A new SupersAttribute object with the result of the subtraction.
        """
        if isinstance(other, SupersAttribute):
            new_data = self.data - other.data
        else:
            new_data = self.data - other
        return SupersAttribute(
            name=self.name,
            data=new_data,
            units=self.units,
            metadata=self.metadata,
        )

    def __mul__(self, other):
        """
        This method overloads the * operator.
        It performs element-wise multiplication of the data attribute with another object.

        Parameters
        ----------
        other : object
            The object to be multiplied with the data attribute.

        Returns
        -------
        SupersAttribute
            A new SupersAttribute object with the result of the multiplication.
        """
        if isinstance(other, SupersAttribute):
            new_data = self.data * other.data
        else:
            new_data = self.data * other
        return SupersAttribute(
            name=self.name,
            data=new_data,
            units=self.units,
            metadata=self.metadata,
        )

    def __truediv__(self, other):
        """
        This method overloads the / operator.
        It performs element-wise division of the data attribute by another object.

        Parameters
        ----------
        other : object
            The object to divide the data attribute by.

        Returns
        -------
        SupersAttribute
            A new SupersAttribute object with the result of the division.
        """
        if isinstance(other, SupersAttribute):
            new_data = self.data / other.data
        else:
            new_data = self.data / other
        return SupersAttribute(
            name=self.name,
            data=new_data,
            units=self.units,
            metadata=self.metadata,
        )

    def __pow__(self, other):
        """
        This method overloads the ** operator.
        It performs element-wise exponentiation of the data attribute by another object.

        Parameters
        ----------
        other : object
            The object to exponentiate the data attribute by.

        Returns
        -------
        SupersAttribute
            A new SupersAttribute object with the result of the exponentiation.
        """
        if isinstance(other, SupersAttribute):
            new_data = self.data**other.data
        else:
            new_data = self.data**other
        return SupersAttribute(
            name=self.name,
            data=new_data,
            units=self.units,
            metadata=self.metadata,
        )

    def attribute_to_indexer(
        self: "SupersAttribute", new_name: Union[str, None] = None
    ) -> "SupersIndexer":
        """
        This function converts an attribute to an indexer.
        The attribute is converted to an indexer by creating a SupersIndexer object.

        Parameters
        ----------
        attribute : SupersAttribute
            The attribute to be converted to an indexer.
        new_name : str, optional
            The new name of the indexer.
            Default is None. If None, the original name is used.

        Returns
        -------
        SupersIndexer
            The indexer created from the attribute.
        """

        if new_name is None:
            new_name = self.name

        return SupersIndexer(
            name=new_name,
            data=self.data,
            units=self.units,
            metadata=self.metadata,
        )

    def attribute_to_indexer_unique(
        self: "SupersAttribute", new_name: Union[str, None] = None
    ) -> "SupersIndexer":
        """
        This function converts an attribute to an indexer.
        The attribute is converted to an indexer by creating a SupersIndexer object.

        Parameters
        ----------
        attribute : SupersAttribute
            The attribute to be converted to an indexer.
        new_name : str, optional
            The new name of the indexer.
            Default is None. If None, the original name is used.

        Returns
        -------
        SupersIndexer
            The indexer created from the attribute.
        """
        if new_name is None:
            new_name = self.name

        return SupersIndexerUnique(
            name=new_name,
            data=self.data,
            units=self.units,
            metadata=self.metadata,
        )

    def attribute_to_indexer_binned(
        self: "SupersAttribute",
        bins: np.ndarray,
        new_name: Union[str, None] = None,
        right: bool = False,
    ) -> "SupersIndexerBinned":
        """
        This function converts an attribute to a binned indexer.
        The attribute is converted to a binned indexer by creating a SupersIndexerBinned object.

        Parameters
        ----------
        attribute : SupersAttribute
            The attribute to be converted to a binned indexer.
        bins : np.ndarray
            The bin edges of the binned indexer.
        new_name : str, optional
            The new name of the binned indexer.
            Default is None. If None, the original name is used.
        right : bool, optional
            As from the numpy documentation:
            Indicating whether the intervals include the right or the left bin edge.
            Default behavior is (right==False) indicating that the interval does not include the right edge.

        Returns
        -------
        SupersIndexerBinned
            The binned indexer created from the attribute.
        """

        if new_name is None:
            new_name = self.name

        return SupersIndexerBinned(
            name=new_name,
            data=self.data,
            bins=bins,
            right=right,
            units=self.units,
            metadata=self.metadata,
        )

    def bin_attribute_by_counts(
        self: "SupersAttribute", counts: ak.Array
    ) -> "SupersAttribute":
        """
        This function bins an attribute by counts.
        The attribute is binned by counts by creating a new attribute with the binned data.

        Parameters
        ----------
        self : SupersAttribute
            The attribute to be binned by counts.
        counts : ak.Array
            The counts to bin the attribute by.
        """

        ndim = self.data.ndim

        if ndim == 1:
            binned_data = sdtracing.binning_by_1D_counts(
                data=self.data,
                counts=counts,
            )
        elif ndim == 2:
            binned_data = sdtracing.binning_by_2D_counts(
                data=self.data,
                counts=counts,
            )
        elif ndim == 3:
            binned_data = sdtracing.binning_by_3D_counts(
                data=self.data,
                counts=counts,
            )

        self.set_data(binned_data)

    def sort_by(self, sort_array: ak.Array) -> "SupersAttribute":
        """
        This function sorts the attribute by a sort array.
        The attribute is sorted by a sort array by creating a new attribute with the sorted data.

        Parameters
        ----------
        self : SupersAttribute
            The attribute to be sorted.
        sort_array : ak.Array
            The array to sort the attribute by.
        """

        # sort the data
        sorted_data = self.data[sort_array]

        self.set_data(sorted_data)


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

        self.make_coord()
        self.make_digitized_data()

    def set_digitized_data(self, digitized_data: Union[ak.Array, np.ndarray]):
        """
        This function sets the data of the indexer.
        The data is stored in the attribute data.

        Parameters
        ----------
        data : ak.Array or np.ndarray
            The data of the indexer.
        """

        # make sure the data has the same shape as the original data
        sdtracing.assert_same_shape(self.data, digitized_data)

        if isinstance(digitized_data, np.ndarray):
            self.digitized_data = ak.Array(digitized_data)
        elif isinstance(digitized_data, ak.Array):
            self.digitized_data = digitized_data
        else:
            raise ValueError("Data must be an ak.Array or np.ndarray")

    def get_digitized_data(self):
        """
        This function returns the digitized data of the attribute.

        Returns
        -------
        ak.Array
            The digitized data of the attribute.
        """
        return self.digitized_data

    def make_digitized_data(self):
        """
        This function sets the digitized data of the indexer.
        The digitized data is stored in the attribute digitized_data.
        In this class, the digitzed data is the same as the data.
        So the indexer should be integer values.
        """

        # digitize the data
        self.digitized_data = self.data

    def make_coord(self):
        """
        This function sets the coord data of the indexer.
        The coord data is stored in the attribute coord.
        In this class, the coord data is the same as the unique values of the data.
        So the indexer should be integer values.
        """

        # digitize the data
        self.set_coord(coord=np.unique(sdtracing.ak_flatten_full(self.data)))

    def set_coord(self, coord: ak.Array):
        """
        This function sets the coord data of the indexer.
        The coord data is stored in the attribute coord.

        Parameters
        ----------
        coord : np.ndarray
            The coord data of the indexer.
        """

        self.coord = coord

    def get_data(self):
        """
        This function returns the data of the attribute.

        Returns
        -------
        ak.Array
            The data of the attribute.
        """
        return self.digitized_data

    def indexer_to_indexer_binned(
        self, bins: np.ndarray, right: bool = False
    ) -> "SupersIndexerBinned":
        """
        This function converts an indexer to a binned indexer.
        The indexer is converted to a binned indexer by creating a SupersIndexerBinned object.

        Parameters
        ----------
        attribute : SupersAttribute
            The attribute to be converted to a binned indexer.
        bins : np.ndarray
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
            name=self.name,
            data=self.data,
            bins=bins,
            right=right,
            units=self.units,
            metadata=self.metadata,
        )

    def indexer_to_attribute(self: "SupersIndexer") -> "SupersAttribute":
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
            name=self.name,
            data=self.data,
            units=self.units,
            metadata=self.metadata,
        )

    def __str__(self):
        """
        This function returns a string representation of the attribute.
        The string representation contains the name of the attribute and the units.

        Returns
        -------
        str
            The string representation of the attribute.
        """
        return f"{self.name} ({self.units})\ncoord: {self.coord}\n{self.data.typestr}\n{self.digitized_data.typestr}"

    def bin_attribute_by_counts(
        self: "SupersAttribute", counts: ak.Array
    ) -> "SupersAttribute":
        """
        This function bins an attribute by counts.
        The attribute is binned by counts by creating a new attribute with the binned data.

        Parameters
        ----------
        self : SupersAttribute
            The attribute to be binned by counts.
        counts : ak.Array
            The counts to bin the attribute by.
        """

        ndim = self.data.ndim

        if ndim == 1:
            binned_data = sdtracing.binning_by_1D_counts(
                data=self.data,
                counts=counts,
            )
            binned_data_digitized = sdtracing.binning_by_1D_counts(
                data=self.digitized_data,
                counts=counts,
            )
        elif ndim == 2:
            binned_data = sdtracing.binning_by_2D_counts(
                data=self.data,
                counts=counts,
            )
            binned_data_digitized = sdtracing.binning_by_2D_counts(
                data=self.digitized_data,
                counts=counts,
            )
        elif ndim == 3:
            binned_data = sdtracing.binning_by_3D_counts(
                data=self.data,
                counts=counts,
            )
            binned_data_digitized = sdtracing.binning_by_3D_counts(
                data=self.digitized_data,
                counts=counts,
            )

        self.set_data(binned_data)
        self.set_digitized_data(binned_data_digitized)

    def sort_by(self, sort_array: ak.Array) -> "SupersAttribute":
        """
        This function sorts the attribute by a sort array.
        The attribute is sorted by a sort array by creating a new attribute with the sorted data.

        Parameters
        ----------
        self : SupersAttribute
            The attribute to be sorted.
        sort_array : ak.Array
            The array to sort the attribute by.
        """

        # sort the data
        sorted_data = self.data[sort_array]
        sorted_digitized_data = self.digitized_data[sort_array]

        self.set_data(sorted_data)
        self.set_digitized_data(sorted_digitized_data)


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
    - bins (np.ndarray): The bin edges of the binned indexer.
    - bin_centers (np.ndarray): The bin centers of the binned indexer.
    """

    def __init__(
        self,
        name: str,
        data: Union[xr.Dataset, np.ndarray, ak.Array],
        bins: np.ndarray,
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
        bins : np.ndarray
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

        self.set_bins(bins=bins, right=right)
        super().__init__(name=name, data=data, units=units, metadata=metadata)

        # make sure, that the digitized data is correct
        # with the following function, we can make sure that the digitized data only contains necessary values and the coords fit to it.
        self.modify_digitized_data()

    def set_bins(self, bins: np.ndarray, right: bool = False):
        """
        This function sets the bins of the binned indexer.
        The bins are stored in the attribute bins.
        The bin centers are stored in the attribute bin_centers.

        Parameters
        ----------
        bins : np.ndarray
            The bin edges of the binned indexer.
        """

        # set the bins
        self.bins = bins

        # create tuples of the bin edges
        bin_edges = [-np.inf] + list(bins) + [np.inf]
        bin_edges = [
            (bin_edges[i], bin_edges[i + 1]) for i in range(len(bin_edges) - 1)
        ]
        # create the bin centers as the mean of the bin edges
        bin_centers = [np.mean(bin_edges[i]) for i in range(len(bin_edges) - 1)]

        self.bin_edges = bin_edges
        self.bin_centers = np.array(bin_centers)
        self.right = right

    def modify_digitized_data(self):
        """
        If the digitized data does not include any exceeding values from the binnning process, we can neglect these indexes.
        """

        min_digitized = ak.min(self.get_digitized_data())
        diff = min_digitized - 0
        # for the lower bound, we remove the first bin center
        if diff > 0:
            # remove the values in the digitized data which are not in use
            self.set_digitized_data(digitized_data=self.get_digitized_data() - diff)
            # remove the exceeding bin centers and edges
            self.bin_centers = self.bin_centers[diff:]
            self.bin_edges = self.bin_edges[diff:]

        # for the higher bound, we simply remove the bin centers which aren't in use
        max_digitized = ak.max(self.get_digitized_data())
        diff = len(self.bin_centers) - max_digitized

        if diff > 0:
            # remove the exceeding bin centers and edges
            self.bin_centers = self.bin_centers[: -diff + 1]
            self.bin_edges = self.bin_edges[: -diff + 1]

        # make sure the coords fit
        self.make_coord()

    def make_coord(self):
        """
        The coord data of the indexer is the bin centers.
        """

        self.set_coord(coord=self.bin_centers)

    def make_digitized_data(self):
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
            self.digitized_data = np.digitize(
                x=self.data, bins=self.bins, right=self.right
            )
        elif self.data.ndim == 2:
            self.digitized_data = sdtracing.ak_digitize_2D(
                x=self.data, bins=self.bins, right=self.right
            )
        elif self.data.ndim == 3:
            self.digitized_data = sdtracing.ak_digitize_3D(
                x=self.data, bins=self.bins, right=self.right
            )
        else:
            raise NotImplementedError(
                "Only 1D, 2D and 3D arrays are supported for digitization till now."
            )


class SupersIndexerUnique(SupersIndexer):
    """
    This class is used to store a binned indexers attribute of the superdroplets dataset.
    It is a subclass of the SupersAttribute class.
    The coords are the unique values of the data.

    It has the following attributes:
    - name (str): The name of the binned indexer.
    - units (str): The units of the binned indexer.
    - data (ak.Array): The data of the binned indexer.
    - digitized_data (ak.Array): The data of the binned indexer as digitized values.
    - metadata (dict): The metadata of the binned indexer.
    - coord (np.ndarray): The unique values of the binned indexer / the values which can be related to the digitized data.
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
        bins : np.ndarray
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

        self.set_data(data)

        super().__init__(name=name, data=data, units=units, metadata=metadata)

        # self.set_digitized_data()

    def make_coord(self):
        """
        This function sets the coord data of the indexer.
        The coord data is stored in the attribute coord.
        In this class, the coord data is the same as the data.
        So the indexer should be integer values.
        """

        # digitize the data
        self.set_coord(coord=np.unique(sdtracing.ak_flatten_full(self.data)))

    def make_digitized_data(self):
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
            digitized_data = np.digitize(x=self.data, bins=self.coord, right=False)
        elif self.data.ndim == 2:
            digitized_data = sdtracing.ak_digitize_2D(
                x=self.data, bins=self.coord, right=False
            )
        elif self.data.ndim == 3:
            digitized_data = sdtracing.ak_digitize_3D(
                x=self.data, bins=self.coord, right=False
            )
        else:
            raise NotImplementedError(
                "Only 1D, 2D and 3D arrays are supported for digitization till now."
            )

        # as the unique values are used as bins, there should be no exceeding values
        # thus, we can remove the 0 bin and the last bin
        digitized_data = digitized_data - 1

        self.set_digitized_data(digitized_data=digitized_data)


class SupersDataNew(SuperdropProperties):
    """
    The class is used to store the superdroplets dataset.
    It is a subclass of the SuperdropProperties class.

    It contains the following attributes:
    - ds (xr.Dataset): The superdroplets dataset.
    - raggedcount (np.ndarray): The ragged count variable.
    - attributes (dict): The attributes of the superdroplets dataset.
        - The keys are the names of the attributes.
        - The values are SupersAttribute objects.

    Note
    ----
    The attributes are stored in the attribute attributes.
    A list of the attribute names is provided in the attribute ``attribute_names``.
    It contains the following variables:
    - sdId
    - sdgbxindex
    - xi
    - radius
    - msol
    - coord3
    - coord1
    - coord2

    The variable time is stored in the attribute time and is created by a seperate constructor function set_time.

    """

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

        Parameters
        ----------
        dataset : os.PathLike or xr.Dataset
            The path to the superdroplets dataset or the dataset itself.
        consts : dict
            The constants of the superdroplets dataset.
        """

        SuperdropProperties.__init__(self, consts=consts)

        self.ds = self.tryopen_dataset(dataset)
        self.raggedcount = self.ds["raggedcount"].values  # ragged count variable

        self.attributes = dict()
        self.set_attributes_from_ds(attribute_names=self.attribute_names)
        self.set_indexes(indexes=[])

        # set mass attribute
        self.set_attribute(
            SupersAttribute(
                name="mass",
                units="kg",
                data=1e-3 * self.mass(self.get_data("radius"), self.get_data("msol")),
            )
        )

        # set represented mass attribute
        mass_represented = self["mass"] * self["xi"]
        mass_represented.set_name("mass_represented")
        self.set_attribute(mass_represented)

    def tryopen_dataset(self, dataset: Union[os.PathLike, xr.Dataset]) -> xr.Dataset:
        if isinstance(dataset, str):
            print("supers dataset: ", dataset)
            return xr.open_dataset(dataset, engine="zarr", consolidated=False)
        elif isinstance(dataset, xr.Dataset):
            return dataset
        else:
            raise ValueError("dataset must be a path or a xarray dataset")

    def set_attributes_from_ds(self, attribute_names):
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
            self.set_attribute_from_ds(name=name)

        # set time attribute
        self.set_time()

    def set_attribute_from_ds(self, name: str):
        """
        This function sets an attribute of the superdroplets dataset.
        The attribute is stored in the attribute attributes.

        Parameters
        ----------
        name : str
            The name of the attribute to be stored.
        """

        try:
            self.attributes[name] = SupersAttribute(name=name, data=self.ds[name])
        except KeyError:
            print(f"Attribute {name} not found in dataset")

    def set_attribute(
        self, attribute: Union[SupersAttribute, SupersIndexer, SupersIndexerBinned]
    ):
        """
        This function sets an attribute of the superdroplets dataset.
        """
        self.attributes[attribute.name] = attribute

    def set_time(self):
        """
        This function sets the time attribute of the superdroplets dataset.
        The time attribute is stored in the attribute time.
        """

        time_attribute = SupersAttribute(
            name="time",
            data=ak.Array(
                np.repeat(
                    self.ds.time.data,
                    self.raggedcount,
                )
            ),
            units=self.ds.time.units,
            metadata=self.ds.time.attrs,
        )
        self.set_attribute(attribute=time_attribute)

    def set_indexes(self, indexes=Tuple["SupersIndexer", "SupersIndexerBinned"]):
        """
        This function sets the indexes of the superdroplets dataset.
        The indexes are stored in the attribute indexes as a tuple.

        Parameters
        ----------
        indexes : Tuple[SupersIndexer, SupersIndexerBinned]
            The indexes to be stored.

        """

        indexes_list = list()
        indexes_name_list = list()

        for index in indexes:
            if isinstance(index, SupersIndexer) or isinstance(
                index, SupersIndexerBinned
            ):
                indexes_list.append(index)
                indexes_name_list.append(index.name)
            else:
                raise ValueError(
                    "indexes must be a SupersIndexer or a SupersIndexerBinned"
                )

        self.indexes = tuple(indexes_list)
        self.indexes_names = tuple(indexes_name_list)

    def add_index(self, index: Union["SupersIndexer", "SupersIndexerBinned"]):
        """
        This function adds an index to the superdroplets dataset.
        The index is stored in the attribute indexes.

        Parameters
        ----------
        index : SupersIndexer or SupersIndexerBinned
            The index to be added.
        """

        if self.is_index(index.name):
            raise ValueError(f"The index '{index.name}' is already set.")
        if isinstance(index, SupersIndexer) or isinstance(index, SupersIndexerBinned):
            self.indexes = self.indexes + (index,)
            self.indexes_names = self.indexes_names + (index.name,)
        else:
            raise ValueError("index must be a SupersIndexer or a SupersIndexerBinned")

    def is_index(self, index: str) -> bool:
        """
        This function checks if an index is already set in the superdroplets dataset.

        Parameters
        ----------
        index : str
            The name of the index to be checked.

        Returns
        -------
        bool
            True if the index is already set, False otherwise.
        """

        return index in self.indexes_names

    def get_data(self, key):
        try:
            return self.attributes[key].get_data()
        except KeyError:
            err = "no known return provided for " + key + " key"
            raise ValueError(err)

    def __getitem__(self, key):
        try:
            return self.attributes[key]
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
        result += "Attributes:\n--------------\n"
        for name, attribute in self.attributes.items():
            result += str(attribute) + "\n"

        result += "\nIndexes:\n--------------\n"
        for index in self.indexes:
            result += index.name + "\n"
            result += str(index.coord) + "\n"

        return result

    def index_by_attribute(self, attribute: Union["SupersAttribute"]):
        """
        This function indexes all attributes of the superdroplets dataset by one attribute.
        The attribute used for indexing is provided by the index_name.
        It will be converted to an indexer

        Note
        ----
        For the binning of the data, the ``digitzied_data`` value of the indexer is used.
        Depending on the dimensionality of the data, the digitized data is done using the numpy digitize function or the sdtracing.ak_digitize_2D or sdtracing.ak_digitize_3D function.


        Parameters
        ----------
        index : SupersAttribute
            The attribute to be used as a new index.

        """

        index = SupersIndexer.attribute_to_indexer(attribute)
        self.index_by_indexer(index=index)

    def index_by_indexer(self, index: Union["SupersAttribute"]):
        """
        This function indexes all attributes of the superdroplets dataset by one attribute.
        The attribute used for indexing is provided by the index_name.
        It will be converted to an indexer

        Note
        ----
        For the binning of the data, the ``digitzied_data`` value of the indexer is used.
        Depending on the dimensionality of the data, the digitized data is done using the numpy digitize function or the sdtracing.ak_digitize_2D or sdtracing.ak_digitize_3D function.


        Parameters
        ----------
        indexer : SupersAttribute
            The attribute to be used as a new index.

        """

        # validate, if the indexer is already used
        if self.is_index(index.name):
            raise ValueError(
                f"The index '{index.name}' is already set.\n Another indexing by it is not possible!"
            )

        ndim = index.data.ndim

        if ndim == 1:
            counts = sdtracing.create_counts_1D(index.digitized_data)
        elif ndim == 2:
            counts = sdtracing.create_counts_2D(index.digitized_data)
        elif ndim == 3:
            counts = sdtracing.create_counts_3D(index.digitized_data)
        else:
            raise ValueError("Indexer must be 1D, 2D or 3D")

        # sort all attributes by the index
        argsort = ak.argsort(index.digitized_data)
        for attr_name in self.attributes:
            attr = self[attr_name]
            attr.sort_by(sort_array=argsort)

        # add the index to the attributes
        self.set_attribute(attribute=index)

        # bin all indexer attributes
        for attr_name in self.attributes:
            attr = self[attr_name]

            attr.bin_attribute_by_counts(counts=counts)

        self.add_index(index=index)

    def attribute_to_DataArray(self, attribute_name):
        """
        This function converts an attribute to a DataArray.
        The attribute is converted to a DataArray by creating a DataArray object.

        Parameters
        ----------
        attribute : SupersAttribute
            The attribute to be converted to a DataArray.

        Returns
        -------
        xr.DataArray
            The DataArray created from the attribute.
        """
        attribute = self[attribute_name]

        data = attribute.data
        shape = sdtracing.get_awkward_shape(data)

        number_of_variable_axis = np.sum(np.isnan(shape))

        if number_of_variable_axis == 0:
            last_dim_variable = False
        elif number_of_variable_axis == 1:
            try:
                sdtracing.assert_only_last_axis_variable(data)
                last_dim_variable = True
            except ValueError:
                raise ValueError("Only the last axis can be variable")
        else:
            raise ValueError("Only one variable axis is allowed at the last position")

        coords = dict()
        dims = list()
        for index in self.indexes:
            dims.append(index.name)
            coord = index.coord
            if isinstance(coord, np.ndarray):
                coords[index.name] = coord
            elif isinstance(coord, ak.Array):
                coords[index.name] = ak.to_numpy(coord)
            else:
                raise ValueError("coord must be a np.ndarray or ak.Array")

        if last_dim_variable is False:
            data = ak.to_numpy(data)
        else:
            data = sdtracing.ak_ragged_to_padded(data)

        data_shape = data.shape

        while len(data_shape) > len(dims):
            new_dim_idx = len(dims)
            new_dim_len = data_shape[new_dim_idx]
            new_dim_name = f"ragged_dimension_{new_dim_idx}"
            dims.append(new_dim_name)
            coords[new_dim_name] = np.arange(new_dim_len)

        return xr.DataArray(
            data=data,
            dims=dims,
            coords=coords,
            name=attribute.name,
            attrs=attribute.metadata,
        )

    def attribute_to_DataArray_reduction(
        self, attribute_name: str, reduction_func: Callable, kwargs: dict = {}
    ):
        """
        This function converts an attribute to a DataArray.
        The attribute is converted to a DataArray by creating a DataArray object.
        If the last dimension is variable, the data will be reduced in dimensions by the reduction function.
        The ``reduction_func`` must be a function which can be applied to the data.
        Usually this should be an awkward array reduction function like ak.sum, ak.mean, ak.min, ak.max, ...


        Parameters
        ----------
        attribute : SupersAttribute
            The attribute to be converted to a DataArray.
        reduction_func : Callable
            The function to reduce the data in the last dimension.
        kwargs : dict, optional
            The keyword arguments for the reduction function.
            Default is an empty dictionary.

        Returns
        -------
        xr.DataArray
            The DataArray created from the attribute.
            And the data is reduced in the last dimension.
        """
        attribute = self[attribute_name]

        data = attribute.data
        shape = sdtracing.get_awkward_shape(data)

        number_of_variable_axis = np.sum(np.isnan(shape))

        if number_of_variable_axis == 0:
            last_dim_variable = False
        elif number_of_variable_axis == 1:
            try:
                sdtracing.assert_only_last_axis_variable(data)
                last_dim_variable = True
            except ValueError:
                raise ValueError("Only the last axis can be variable")
        else:
            raise ValueError("Only one variable axis is allowed at the last position")

        coords = dict()
        dims = list()
        for index in self.indexes:
            dims.append(index.name)
            coord = index.coord
            if isinstance(coord, np.ndarray):
                coords[index.name] = coord
            elif isinstance(coord, ak.Array):
                coords[index.name] = ak.to_numpy(coord)
            else:
                raise ValueError("coord must be a np.ndarray or ak.Array")

        if last_dim_variable is False:
            data = ak.to_numpy(data)
        else:
            data = reduction_func(data, axis=-1, **kwargs)
            shape = sdtracing.get_awkward_shape(data)
            if any(np.isnan(shape)):
                raise ValueError(
                    "The data is not reduced in a form that no variable axis is left"
                )
            else:
                data = ak.to_numpy(data)

        data_shape = data.shape

        while len(data_shape) > len(dims):
            new_dim_idx = len(dims)
            new_dim_len = data_shape[new_dim_idx]
            new_dim_name = f"ragged_dimension_{new_dim_idx}"
            dims.append(new_dim_name)
            coords[new_dim_name] = np.arange(new_dim_len)

        attrs = attribute.metadata
        attrs.update({"units": attribute.units})

        return xr.DataArray(
            data=data,
            dims=dims,
            coords=coords,
            name=attribute.name,
            attrs=attribute.metadata,
        )
