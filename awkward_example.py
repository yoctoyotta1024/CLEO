# %%
import awkward as ak
import numpy as np
from pySD.sdmout_src import pysetuptxt, pyzarr
from typing import Union
import pytest

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


def eulerian_from_count(data, counts, G):
    # def eulerian_from_counts(data, counts, G) :
    res = ak.unflatten(data, ak.flatten(counts), axis=1)
    # flatten this array again
    res = ak.unflatten(ak.flatten(res), G)
    # then unflatten with the correct length of the number of gridboxes
    res = ak.unflatten(ak.flatten(res), G)
    return res


def create_counts_fast(data: IntegerArray) -> ak.highlevel.Array:
    """
    This function creates a 2D array with the counts of unique values in the input array.
    It identifies the unique values in each subarray along axis 1.

    Parameters
    ----------
    data : IntegerArray
        The input array containing the data for which counts need to be calculated.
        Only ``int`` types are allowed!
        If an ``np.ndarray`` is used, the dimensions should be T x N
        If an ``ak.highlevel.Array`` is used, the dimensions should be T x var

    Returns
    -------
    ak.highlevel.Array
        A 2D array with the counts of unique values in the input array.
    """

    # check for int type:

    assert (
        data.ndim == 2
    ), f"data needs to be 2 dimensional but is {data.ndim} dimensional"

    if isinstance(data, np.ndarray):
        dtype = data.dtype
    elif isinstance(data, ak.highlevel.Array):
        dtype = data.type.content.content.primitive
    else:
        raise TypeError(
            f"The provided type of 'data': '{type(data)}' is not supported!"
        )

    if not np.issubdtype(dtype, np.integer):
        raise TypeError(
            f"The array must contain 'int' like dtypes but has '{dtype}'. This is not supported"
        )

    M = int(ak.max(data)) + 1
    T = int(ak.num(data, axis=0))
    # modified = data + np.arange(T) * M
    # modified = ak.flatten(modified)

    counts_flat = np.bincount(ak.flatten(data + np.arange(T) * M))
    counts_flat = ak.fill_none(ak.pad_none(counts_flat, M * T, axis=0), 0)
    counts_flat = ak.unflatten(counts_flat, M)
    return counts_flat


def euler_full(
    data: ak.highlevel.Array, indexer: ak.highlevel.Array
) -> ak.highlevel.Array:
    """
    Calculates the Eulerian from the given data and indexer.
    The indexer array needs to have a ``int`` like dtype.
    The output shape of the array will be given by the maximum value M in ``indexer``.
    Output shape : T x M x var

    Note
    ----
    This functions uses np.bincount on the lower levels.
    Thus, at one point an np.ndarray of shape (T * M) will be created and needs to be stored in memory.
    This limits the capability due to memory usage.

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
    N = int(ak.max(indexer)) + 1

    # sort arrays by their
    argsort = ak.argsort(indexer, axis=1)
    indexer_sort = indexer[argsort]
    data_sort = data[argsort]

    counts = create_counts_fast(indexer_sort)

    # The function eulerian_from_count performes this bit of code:
    # ----------
    # # unflatten the data array using the counts array.
    # res = ak.unflatten(data, ak.flatten(counts), axis=1)
    # # flatten this array again
    # res = ak.unflatten(ak.flatten(res), N)
    # # then unflatten with the correct length of the number of gridboxes
    # res = ak.unflatten(ak.flatten(res), N)
    # return res
    # ----------

    return eulerian_from_count(data_sort, counts, N)


# %%

data = ak.Array(
    [
        [10.2, 20, 30],
        [12],
        [40, 50],
        [90],
    ]
)
time_index = ak.Array(
    [
        0,
        1,
        2,
        3,
    ]
)
gridbox = ak.Array(
    [
        [5, 1, 2],
        [1],
        [3, 3],
        [3],
    ]
)

counts = ak.Array(
    [
        [2, 1, 0, 0],
        [0, 0, 0, 0],
        [0, 1, 1, 0],
        [0, 0, 0, 1],
    ]
)

id = ak.Array(
    [
        [3, 1, 2],
        [2],
        [1, 3],
        [7],
    ]
)

setupfile = "/home/m/m301096/CLEO/data/output/raw/no_aerosols_collision_many_5012/clusters_301/eurec4a1d_setup.txt"
dataset = "/home/m/m301096/CLEO/data/output/raw/no_aerosols_collision_many_5012/clusters_301/eurec4a1d_sol.zarr"
# read in constants and intial setup from setup .txt file
config = pysetuptxt.get_config(setupfile, nattrs=3, isprint=False)
consts = pysetuptxt.get_consts(setupfile, isprint=False)
# Create a first simple dataset to have the coordinates for later netcdf creation
sddata = pyzarr.get_supers(dataset, consts)


# data = sddata.radius
# id = sddata.sdId + int(10e5)
# gridbox = sddata.sdgbxindex

euler_full(data, id)


# def calc_memory_array(elements):
#     memory_usage = elements * 8  # 8 bytes per float

#     if memory_usage < 1024 * 1024:
#         memory_usage_mb = memory_usage / (1024 * 1024)
#         return f"Memory usage: {memory_usage_mb:.2f} Mb"
#     elif memory_usage < 1024 * 1024 * 1024:
#         memory_usage_gb = memory_usage / (1024 * 1024 * 1024)
#         return f"Memory usage: {memory_usage_gb:.2f} Gb"
#     else:
#         memory_usage_tb = memory_usage / (1024 * 1024 * 1024 * 1024)
#         return f"Memory usage: {memory_usage_tb:.2f} Tb"
