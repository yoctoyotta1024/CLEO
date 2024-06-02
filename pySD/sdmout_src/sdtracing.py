"""
----- CLEO -----
File: sdtracing.py
Project: sdmout_src
Created Date: Tuesday 24th October 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
Last Modified: Monday 15th April 2024
Modified By: CB
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
Copyright (c) 2023 MPI-M, Clara Bayley
-----
File Description:
functions to extract attribute data for
specifc superdroplets based on the sdIds
e.g. for tracing their trajectories
"""

# %%
import numpy as np
import awkward as ak
import random
import warnings
from typing import Union


def attr_for_superdroplet(sddata, Id, attr):
    """selects attribute from sddata belonging
    to superdroplet with identitiy 'Id'
    at every output time"""

    bools = ak.Array(sddata.sdId == Id)  # True/False id is found in sdId at each time
    attr4Id = sddata[attr][bools]  # attribute where sdId = Id

    num = ak.num(
        attr4Id
    )  # at each time, number of positions where sdId is found (should be 0 or 1)
    if any(num[num != 1]):
        errmsg = (
            "attribute timeseries has times when more"
            + " than one sdId==Id. num should be list of either 1 or 0"
            + " (for Id found in sddata at given time or not)"
        )
        raise ValueError(errmsg)

    attr4Id = ak.where(
        num == 0, ak.Array([[np.nan]]), attr4Id
    )  # replace empty values with np.nan

    return ak.flatten(attr4Id, axis=1)  # remove excess dimension


def attributes_for1superdroplet(sddata, Id, attrs):
    """selects attributes in 'attrs' from sddata
    belonging to superdroplet with identitiy 'Id'
    at every output time"""

    attrs4Id = {}
    for attr in attrs:
        attrs4Id[attr] = attr_for_superdroplet(sddata, Id, attr)

    return attrs4Id


def attribute_for_superdroplets_sample(
    sddata, attr, ndrops2sample=0, minid=0, maxid=0, ids=[]
):
    """returns 2D array with dimensions [time, SD]
    containing attribute data over time for a sample of
    superdroplets. Sample is either for superdroplets with
    specific Ids in 'ids' list, or sample of 'ndrops2sample'
    randomly selected superdrops with Ids in the range
    [minid, maxid]"""

    if np.any(ids):
        sample = ids
    else:  # ids == []
        population = list(range(minid, maxid, 1))
        if ndrops2sample == 0:
            ndrops2sample = maxid
        sample = random.sample(population, ndrops2sample)

    ndrops_attr = []
    for id in sample:
        attr4Id = attr_for_superdroplet(sddata, id, attr)
        ndrops_attr.append(attr4Id)

    return np.asarray(ndrops_attr).T


def attr_at_times(attrdata, time, times2sel):
    """selects attribute (for all superdroplets)
    at times closest to 'times2sel'"""

    inds = []  # list containing indexes of times closest to times2sel
    for t in times2sel:
        inds.append(np.argmin(abs(time - t)))

    return attrdata[inds]


def attributes_at_times(sddata, time, times2sel, attrs2sel):
    """selects attributes at given times from
    sddata (for all superdroplets in sddata)"""

    selected_data = {}  # dict containting selected attributes at selected times

    for attr in attrs2sel:
        selattr_data = attr_at_times(sddata[attr], time, times2sel)
        selected_data[attr] = selattr_data

    return selected_data


def attrs_for_superdroplets_sample(
    sddata, attrs, ndrops2sample=0, minid=0, maxid=0, ids=[]
):
    """returns dictionary of 2D arrays (with dimensions [time, SD])
    for each attribute in 'attrs' list for a sample of
    superdroplets. Sample is either for superdroplets with
    specific Ids in 'ids' list, or sample of 'ndrops2sample'
    randomly selected superdrops with Ids in the range
    [minid, maxid]"""

    if np.any(ids):
        sample = ids
    else:  # ids == []
        population = list(range(minid, maxid, 1))
        if ndrops2sample == 0:
            ndrops2sample = maxid
        sample = random.sample(population, ndrops2sample)

    data = {}
    for a in attrs:
        data[a] = attribute_for_superdroplets_sample(sddata, a, ids=sample)

    return data


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


def get_awkward_shape(a: ak.highlevel.Array):
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


def assert_same_shape(a: ak.highlevel.Array, b: ak.highlevel.Array):
    """
    Assert that the two awkward arrays have the same shape.
    Variable axis lengths are replaced by ``np.nan``.

    Parameters
    ----------
    a : ak.Array
        The first input array.
    b : ak.Array
        The second input array.

    Raises
    ------
    ValueError
        If the shapes of the two arrays do not match.
    """

    shape_a = get_awkward_shape(a)
    shape_b = get_awkward_shape(b)
    if shape_a != shape_b:
        raise ValueError(
            f"The shapes of the two arrays do not match: {shape_a} != {shape_b}"
        )


def assert_only_last_axis_variable(a: ak.highlevel.Array):
    """
    Assert that the awkward array has only variable axis lengths at the last axis.

    Parameters
    ----------
    a : ak.Array
        The input array.

    Raises
    ------
    ValueError
        If the array has variable axis lengths at other axes than the last axis.
    """

    shape = get_awkward_shape(a)
    if any([np.isnan(elem) for elem in shape[:-1]]):
        raise ValueError(
            "The array has variable axis lengths at other axes than the last axis."
        )


def ak_flatten_full(data: ak.Array) -> ak.highlevel.Array:
    """
    This function flattens the input array along all axes.
    Note that no information is returned to unflatten the array again!

    Parameters
    ----------
    data : ak.highlevel.Array
        The input array to be flattened.

    Returns
    -------
    ak.highlevel.Array
        The flattened array.
    """

    while data.ndim > 1:
        data = ak.flatten(data)

    return data


def ak_ragged_to_padded(
    data: ak.Array, masked_array: bool = False
) -> Union[np.ndarray, np.ma.MaskedArray]:
    """
    This function converts a ragged array to a padded array and returns it as a numpy array.
    If the last axis has variable lengths, the array will be padded to the maximum length of the last axis.
    This behaviour can lead to very large arrays, especilly if the array is very sparse.
    Missing values are filled with np.nan values.

    Setting the ``masked_array`` parameter to ``True`` will return a masked array instead of a regular numpy array.

    Parameters
    ----------
    data : ak.highlevel.Array
        The input ragged array.
    masked_array : bool, optional
        If set to ``True``, a masked array is returned.
        Otherwise a regular numpy array is returned with np.nan values at the missing positions.
        The default is ``False``.

    Returns
    -------
    np.ndarray or np.ma.MaskedArray
        The padded array in numpy format.

    Raises
    ------
    ValueError
        If the input array has variable axis lengths at other axes than the last axis.

    """

    assert_only_last_axis_variable(data)

    pad_max = int(ak.max(ak.num(data, axis=-1)))
    padded_data = ak.pad_none(data, pad_max, axis=-1)

    padded_data_numpy = padded_data.to_numpy().astype(float)
    padded_data_numpy.set_fill_value(np.nan)
    if masked_array is True:
        return padded_data_numpy
    else:
        return padded_data_numpy.filled()


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


def create_counts_1D(a: ak.highlevel.Array) -> ak.highlevel.Array:
    """
    This function creates a 1D array with the counts of unique values in the input array.
    It identifies the unique values in each subarray along axis 1.

    Parameters
    ----------
    a : ak.highlevel.Array
        The input array of shape T, containing the data for which counts need to be calculated.
        Its maximum value ``B`` will be used to identify
        the number of bins as ``M = B + 1``
        Only ``int`` types are allowed!
    Returns
    -------
    ak.highlevel.Array
        An array with the counts of unique values in the input array.
        The output shape will be (M).
    """

    # check for int type:

    assert a.ndim == 1, f"data needs to be 1 dimensional but is {a.ndim} dimensional"

    if isinstance(a, np.ndarray):
        dtype = a.dtype
    elif isinstance(a, ak.highlevel.Array):
        dtype = a.type.content.primitive
    else:
        raise TypeError(f"The provided type of 'data': '{type(a)}' is not supported!")

    if not np.issubdtype(dtype, np.integer):
        raise TypeError(
            f"The array must contain 'int' like dtypes but has '{dtype}'. This is not supported"
        )

    M = int(ak.max(a)) + 1
    # create a bincount array with minimal length of the maximum value in the array
    bcount = np.bincount(a, minlength=M)
    return bcount


def create_counts_2D(a: ak.highlevel.Array, flat: bool = False) -> ak.highlevel.Array:
    """
    This function creates a 2D array with the counts of unique values in the input array.
    It identifies the unique values in each subarray along axis 1.

    Parameters
    ----------
    a : ak.highlevel.Array
        The input array of shape T x var, containing the data for which counts need to be calculated.
        Its maximum value ``B`` will be used to identify
        the number of bins as ``M = B + 1``
        Only ``int`` types are allowed!
    flat : bool, optional
        If False, the output array will have the shape (T, M).
        If True, the output array is flattened top be a 1D array of length T*M.
        Default is False.
    Returns
    -------
    ak.highlevel.Array
        An array with the counts of unique values in the input array.
        The output shape will be (T, M) if flat is False.
        The output shape will be (T*M) if flat is True.

    """

    # check for int type:

    assert a.ndim == 2, f"data needs to be 2 dimensional but is {a.ndim} dimensional"

    if isinstance(a, np.ndarray):
        dtype = a.dtype
    elif isinstance(a, ak.highlevel.Array):
        dtype = a.type.content.content.primitive
    else:
        raise TypeError(f"The provided type of 'data': '{type(a)}' is not supported!")

    if not np.issubdtype(dtype, np.integer):
        raise TypeError(
            f"The array must contain 'int' like dtypes but has '{dtype}'. This is not supported"
        )

    T = int(ak.num(a, axis=0))
    M = int(ak.max(a)) + 1

    # modify the arra to have unique values
    # for each combination of i and j along dim0, dim1
    modified = a + np.arange(T) * M
    modified = ak.flatten(modified)

    bcount = np.bincount(modified, minlength=T * M)

    if flat is False:
        bcount = ak.unflatten(bcount, M)
        return bcount
    else:
        return bcount


def create_counts_3D(a: ak.highlevel.Array, flat: bool = False) -> ak.highlevel.Array:
    """
    This function creates a 3D array with the counts of unique values in the input array.
    For instance, if the ``a`` array has the shape (T, N, var) and its highest values is M, the output array will have the shape (T, N, M+1),

    Parameters
    ----------
    a : ak.highlevel.Array
        The input array of shape T x N x var, containing the data for which counts need to be calculated.
        Its maximum value ``B`` will be used to identify
        the number of bins as ``M = B + 1``
        Only ``int`` types are allowed!
    flat : bool, optional
        If False, the output array will have the shape (T, N, M).
        If True, the output array is flattened top be a 1D array of length T*N*M.
        Default is False.

    Returns
    -------
    ak.highlevel.Array
        An array with the counts of values given in the indexer array.
        The output shape will be (T, N, M) if flat is False.
        The output shape will be (T*N*M) if flat is True.

    Notes
    -----
    To make sure the indexer has unique value for each combination
    of i and j along axis 0 and axis 1, the input array ``a`` is modified, by adding a
    2D array to it. This 2D array is created by np.arange(0, x * y * z, z).reshape(x, y).

    Example with
    T = 4, N = 2, B = 4, M = 5
    x = 4, y = 2, z = 5
    >>> a = ak.Array([
            [[0, 1], [1]],
            [[0],    [1, 4]],
            [[1, 3], [1]],
            [[0],    [3, 4]]
        ])
    >>> add_array = np.arange(0, x * y * z, z).reshape(x, y)
    >>> add_array
    ... [[ 0,  5],
    ...  [10, 15],
    ...  [20, 25],
    ...  [30, 35]]
    >>> modified = a + add_array
    >>> modified
    ... [
    ...  [[ 0,  1],  [ 6]],
    ...  [[10],      [16, 19]],
    ...  [[21, 23],  [26]],
    ...  [[30],      [38, 39]],
    ... ]

    """

    assert a.ndim == 3, f"data needs to be 3 dimensional but is {a.ndim} dimensional"

    if isinstance(a, np.ndarray):
        dtype = a.dtype
    elif isinstance(a, ak.highlevel.Array):
        dtype = a.type.content.content.content.primitive
    else:
        raise TypeError(f"The provided type of 'data': '{type(a)}' is not supported!")

    if not np.issubdtype(dtype, np.integer):
        raise TypeError(
            f"The array must contain 'int' like dtypes but has '{dtype}'. This is not supported"
        )

    T, N, _ = get_awkward_shape(a)
    M = int(ak.max(a)) + 1

    # to make sure the indexer has unique value for each combination
    # of i and j along dim0, dim1, a 2D array is created.
    # Example:
    # T = 4, N = 2, M = 4
    # x = 4, y = 3, z = 5
    # [[ 0,  5],
    #  [10, 15],
    #  [20, 25],
    #  [30, 35]]

    add_array = np.arange(0, T * N * M, M).reshape(T, N)

    modified = a + add_array
    modified = ak_flatten_full(modified)

    bcount = np.bincount(modified, minlength=T * N * M)

    if flat is False:
        bcount = ak.unflatten(bcount, M)
        bcount = ak.unflatten(bcount, N)
        return bcount
    else:
        return bcount


def binning_by_1D_indexer(
    data: ak.highlevel.Array,
    indexer: ak.highlevel.Array,
) -> ak.highlevel.Array:
    """
    Calculates the Eulerian 2D array from the given 1D data and 1D indexer arrays.
    The indexer array needs to have a ``int`` like dtype.
    The output shape of the array will be given by the maximum value M in ``indexer``.
    Output shape : T x M x var

    Note
    ----
    - This functions uses np.bincount on the lower levels.
    - Thus, at one point an np.ndarray of shape (M) will be created and needs to be stored in memory.
    - This limits the capability due to memory usage.
    - A good practice is, to extract the N unique values of the indexer array. With this, create a new indexer with values from 0 to N. You can use the numpy.digitize function for this.

    >>> unique_index = np.unique(indexer)
    >>> unique_indexer = np.digitize(indexer, unique_index)
    >>> binning_by_1D_indexer(data, unique_indexer)

    For a lagrangian tracking for non sparse outputs, it is prefered to use the lagrangian function.

    Parameters
    ----------
    data (ak.highlevel.Array):
        The input data array (T).
    indexer (ak.highlevel.Array):
        The indexer array (T) of dtype int.
        With maximum value M.

    Returns
    -------
    ak.highlevel.Array:
        The calculated Eulerian array of shape (T x M x var)
    """

    assert_same_shape(data, indexer)

    # sort arrays by their
    argsort = ak.argsort(indexer, axis=0)
    indexer_sort = indexer[argsort]
    data_sort = data[argsort]

    # use np.bincount to get the counts of the unique values
    counts = create_counts_1D(indexer_sort)
    # unflatten the data array using the bcount array.
    result = ak.unflatten(data_sort, counts, axis=0)
    return result


def binning_by_2D_indexer(
    data: ak.highlevel.Array, indexer: ak.highlevel.Array
) -> ak.highlevel.Array:
    """
    Calculates the Eulerian 3D array from the given 2D data and 2D indexer arrays.
    The indexer array needs to have a ``int`` like dtype.
    The output shape of the array will be given by the maximum value M in ``indexer``.
    Output shape : T x M x var

    Note
    ----
    - This functions uses np.bincount on the lower levels.
    - Thus, at one point an np.ndarray of shape (T * M) will be created and needs to be stored in memory.
    - This limits the capability due to memory usage.
    - A good practice is, to extract the N unique values of the indexer array. With this, create a new indexer with values from 0 to N. You can use the numpy.digitize function for this.

    >>> unique_index = np.unique(ak.flatten(indexer))
    >>> unique_indexer = ak_digitize_2d(indexer, unique_index)
    >>> binning_by_2D_indexer(data, unique_indexer)

    For a lagrangian tracking for non sparse outputs, it is prefered to use the lagrangian function.

    Parameters
    ----------
    data (ak.highlevel.Array):
        The input data array (T x var).
    indexer (ak.highlevel.Array):
        The indexer array (T x var) of dtype int.
        With maximum value M.

    Returns
    -------
    ak.highlevel.Array:
        The calculated Eulerian array of shape (T x M x var)
    """
    # indexer_unique = np.unique(ak.flatten(indexer))
    assert_same_shape(data, indexer)
    assert_only_last_axis_variable(data)
    assert_only_last_axis_variable(indexer)

    N = int(ak.max(indexer)) + 1

    # sort arrays by their
    argsort = ak.argsort(indexer, axis=1)
    indexer_sort = indexer[argsort]
    data_sort = data[argsort]

    counts = create_counts_2D(indexer_sort)

    # The function eulerian_from_count performes this bit of code:
    # ----------
    # # unflatten the data array using the counts array.
    result = ak.unflatten(data_sort, ak.flatten(counts), axis=1)
    # flatten this array again
    result = ak.unflatten(ak.flatten(result), N)

    return result


def binning_by_3D_indexer(
    data: ak.highlevel.Array, indexer: ak.highlevel.Array
) -> ak.highlevel.Array:
    """
    This function bins the data array according to the indexer arrays values, which should be integers.
    The resulting array will have the shape (T, N, B, var), where B is the number of unique values in the indexer array.

    Note
    ----
    - The input arrays must have the same shape.
    - The indexer array must have integer values.
    - The data array can ONLY have variable axis lengths at the last axis!

    Parameters
    ----------
    data : ak.Array
        The input data array with shape (T, N, var).
    indexer : ak.Array
        The indexer array with shape (T, N, var).
        And B unique values.

    Returns
    -------
    result : ak.highlevel.Array
        The binned data array with shape (T, N, B, var).
    """

    assert_same_shape(data, indexer)
    assert_only_last_axis_variable(data)
    assert_only_last_axis_variable(indexer)

    x, y, _ = get_awkward_shape(indexer)
    z = int(ak.max(indexer)) + 1

    args = ak.argsort(indexer, axis=2)
    data = data[args]
    indexer = indexer[args]

    # flatten the data array
    data_flat = ak.flatten(data, axis=1)
    data_flat = ak.flatten(data_flat, axis=1)

    # to make sure the indexer has unique value for each combination od i and j along dim0, dim1, a 2D array is created.
    # Example:
    # T = 4, N = 3, B = 4
    # x = 4, y = 3, z = 5
    # [[ 0,  5, 10],
    #  [15, 20, 25],
    #  [30, 35, 40],
    #  [45, 50, 55]])

    counts = create_counts_3D(indexer, flat=True)

    result = ak.unflatten(data_flat, counts)
    result = ak.unflatten(result, z)
    result = ak.unflatten(result, y)

    return result


def binning_by_1D_counts(
    data: ak.highlevel.Array,
    counts: ak.highlevel.Array,
) -> ak.highlevel.Array:
    """
    Calculates the Eulerian 2D array from the given 1D data and 1D indexer arrays.
    The indexer array needs to have a ``int`` like dtype.
    The output shape of the array will be given by the maximum value M in ``indexer``.
    Output shape : T x M x var

    Note
    ----
    - This functions uses np.bincount on the lower levels.
    - Thus, at one point an np.ndarray of shape (M) will be created and needs to be stored in memory.
    - This limits the capability due to memory usage.
    - A good practice is, to extract the N unique values of the indexer array. With this, create a new indexer with values from 0 to N. You can use the numpy.digitize function for this.

    >>> unique_index = np.unique(indexer)
    >>> unique_indexer = np.digitize(indexer, unique_index)
    >>> binning_by_1D_indexer(data, unique_indexer)

    For a lagrangian tracking for non sparse outputs, it is prefered to use the lagrangian function.

    Parameters
    ----------
    data (ak.highlevel.Array):
        The input data array (T).
    indexer (ak.highlevel.Array):
        The indexer array (T) of dtype int.
        With maximum value M.

    Returns
    -------
    ak.highlevel.Array:
        The calculated Eulerian array of shape (T x M x var)
    """

    if ak.sum(counts) != ak.count(data):
        raise ValueError(
            "The sum of counts does not match the number of elements in the data array"
        )

    # unflatten the data array using the bcount array.
    result = ak.unflatten(data, counts, axis=0)
    return result


def binning_by_2D_counts(
    data: ak.highlevel.Array, counts: ak.highlevel.Array
) -> ak.highlevel.Array:
    """
    Calculates the Eulerian 3D array from the given 2D data and 2D indexer arrays.
    The indexer array needs to have a ``int`` like dtype.
    The output shape of the array will be given by the maximum value M in ``indexer``.
    Output shape : T x M x var

    Note
    ----
    - This functions uses np.bincount on the lower levels.
    - Thus, at one point an np.ndarray of shape (T * M) will be created and needs to be stored in memory.
    - This limits the capability due to memory usage.
    - A good practice is, to extract the N unique values of the indexer array. With this, create a new indexer with values from 0 to N. You can use the numpy.digitize function for this.

    >>> unique_index = np.unique(ak.flatten(indexer))
    >>> unique_indexer = ak_digitize_2d(indexer, unique_index)
    >>> binning_by_2D_indexer(data, unique_indexer)

    For a lagrangian tracking for non sparse outputs, it is prefered to use the lagrangian function.

    Parameters
    ----------
    data (ak.highlevel.Array):
        The input data array (T x var).
    indexer (ak.highlevel.Array):
        The indexer array (T x var) of dtype int.
        With maximum value M.

    Returns
    -------
    ak.highlevel.Array:
        The calculated Eulerian array of shape (T x M x var)
    """
    # indexer_unique = np.unique(ak.flatten(indexer))

    data_shape = get_awkward_shape(data)

    if len(data_shape) != 2:
        raise ValueError("The data array must be 2D")

    assert_only_last_axis_variable(data)

    counts_shape = get_awkward_shape(counts)
    if len(counts_shape) != 2:
        raise ValueError("The counts array must be 2D")

    N = counts_shape[1]

    if ak.sum(counts) != ak.count(ak.flatten(data)):
        raise ValueError(
            "The sum of counts does not match the number of elements in the data array"
        )

    # # unflatten the data array using the counts array.
    counts_flat = ak_flatten_full(counts)

    result = ak.unflatten(data, counts_flat, axis=1)
    # flatten this array again
    result = ak.unflatten(ak.flatten(result), N)

    return result


def binning_by_3D_counts(
    data: ak.highlevel.Array, counts: ak.highlevel.Array
) -> ak.highlevel.Array:
    """
    This function bins the data array according to the indexer arrays values, which should be integers.
    The resulting array will have the shape (T, N, B, var), where B is the number of unique values in the indexer array.

    Note
    ----
    - The input arrays must have the same shape.
    - The indexer array must have integer values.
    - The data array can ONLY have variable axis lengths at the last axis!

    Parameters
    ----------
    data : ak.Array
        The input data array with shape (T, N, var).
    indexer : ak.Array
        The indexer array with shape (T, N, var).
        And B unique values.

    Returns
    -------
    result : ak.highlevel.Array
        The binned data array with shape (T, N, B, var).
    """

    data_shape = get_awkward_shape(data)

    if len(data_shape) != 3:
        raise ValueError("The data array must be 2D")

    assert_only_last_axis_variable(data)

    counts_shape = get_awkward_shape(counts)
    if len(counts_shape) != 3:
        raise ValueError("The counts array must be 3D")

    x, y, _ = get_awkward_shape(data)
    z = counts_shape[2]

    # flatten the data array
    data_flat = ak.flatten(data, axis=1)
    data_flat = ak.flatten(data_flat, axis=1)

    if ak.sum(counts) != ak.count(data_flat):
        raise ValueError(
            "The sum of counts does not match the number of elements in the data array"
        )

    counts_flat = ak_flatten_full(counts)

    result = ak.unflatten(data_flat, counts_flat)
    result = ak.unflatten(result, z)
    result = ak.unflatten(result, y)

    return result
