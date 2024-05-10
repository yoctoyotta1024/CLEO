# %%
import awkward as ak
import numpy as np
from tqdm import tqdm
from pySD.sdmout_src import pysetuptxt, pyzarr
from typing import Union


IntegerArray = Union[np.ndarray, ak.highlevel.Array]

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

# %%


def eulerian_from_count(data, counts, G):
    # def eulerian_from_counts(data, counts, G) :
    res = ak.unflatten(data, ak.flatten(counts), axis=1)
    # flatten this array again
    res = ak.unflatten(ak.flatten(res), G)
    # then unflatten with the correct length of the number of gridboxes
    res = ak.unflatten(ak.flatten(res), G)
    return res


def create_counts_np(data, counts):
    T, G = np.shape(counts)
    # for each timestep, identify the unique values and their counts
    for t in tqdm(range(T)):
        u, c = np.unique(data[t], return_counts=True, axis=0)
        counts[t, u] = c
    return counts


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
    modified = data + np.arange(T) * M
    modified = ak.flatten(modified)

    counts_flat = np.bincount(modified)
    counts_flat = ak.fill_none(ak.pad_none(counts_flat, M * T, axis=0), 0)
    counts_flat = ak.unflatten(counts_flat, M)
    return counts_flat


# %%

time_unique = np.unique(time_index)
gridbox_unique = np.unique(ak.flatten(gridbox))
id_unique = np.unique(ak.flatten(id))
T = int(ak.num(time_unique, axis=0))
G = int(ak.num(gridbox_unique, axis=0))
N = int(ak.max(id_unique, axis=0)) + 1


gridbox_argsort = ak.argsort(gridbox)
gridbox_argsort = ak.argsort(id)
id_sort = id[gridbox_argsort]
gridbox_sort = gridbox[gridbox_argsort]
data_sort = data[gridbox_argsort]
id_sort = id[gridbox_argsort]

# %%


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


euler_full(data, id)
# %%


def cut_unique_values(array, data):
    if len(array) != 0:
        unique, counts = np.unique(array, return_counts=True, axis=0)
        return ak.unflatten(data, counts)
    else:
        return array


# %%
def empty_list_of_lists(length=[10, 20]):
    dim = len(length)
    result = []
    for i in range(length[0]):
        result.append(empty_list_of_lists(length[1:]) if dim > 1 else [])
    return result


def euler(data_sort, gridbox_sort, T, G):
    result = empty_list_of_lists([T, G])
    # result = np.array((T, G), dtype = object)
    for t in range(T):
        data_now = data_sort[t]
        grid_now = gridbox_sort[t]
        if len(data_now) == 0:
            pass
        else:
            grid_split = cut_unique_values(grid_now, grid_now)
            data_split = cut_unique_values(grid_now, data_now)

            grid_split_flat = ak.flatten(grid_split)
            data_split_flat = ak.flatten(data_split)

            for i in range(len(grid_split_flat)):
                result[t][grid_split_flat[i]].append(data_split_flat[i])

    return ak.Array(result)


def euler2(data_sort, gridbox_sort, T, G):
    result = empty_list_of_lists([T, G])
    for t in tqdm(range(T)):
        data_now = data_sort[t]
        grid_now = gridbox_sort[t]
        if len(data_now) == 0:
            pass
        else:
            grid_split = cut_unique_values(grid_now, grid_now)
            data_split = cut_unique_values(grid_now, data_now)
            for i in range(len(grid_split)):
                result[t][grid_split[i][0]] = data_split[i]

    return result


# %%

setupfile = "/home/m/m301096/CLEO/data/output/raw/no_aerosols_collision_many_5012/clusters_301/eurec4a1d_setup.txt"
dataset = "/home/m/m301096/CLEO/data/output/raw/no_aerosols_collision_many_5012/clusters_301/eurec4a1d_sol.zarr"
# read in constants and intial setup from setup .txt file
config = pysetuptxt.get_config(setupfile, nattrs=3, isprint=False)
consts = pysetuptxt.get_consts(setupfile, isprint=False)
# Create a first simple dataset to have the coordinates for later netcdf creation
sddata = pyzarr.get_supers(dataset, consts)


data = sddata.radius
id = sddata.sdId + int(1e5)
gridbox = sddata.sdgbxindex

time_index = sddata.time
time_unique = np.unique(sddata.time)
gridbox_unique = np.unique(ak.flatten(sddata.sdgbxindex))
id_unique = np.unique(ak.flatten(sddata.sdId))

T = int(ak.num(time_unique, axis=0))
G = int(ak.num(gridbox_unique, axis=0))
N = int(ak.max(id_unique)) + 1

gridbox_argsort = ak.argsort(sddata.sdgbxindex)
gridbox_argsort = ak.argsort(sddata.sdId)
gridbox_sort = sddata.sdgbxindex[gridbox_argsort]
id_sort = sddata.sdId[gridbox_argsort]
data_sort = sddata.radius[gridbox_argsort]

# %%
# counts_fast = create_counts_fast(id_sort, T, G)
# eulerian_from_count(data_sort, create_counts_fast(gridbox_sort, T, G), G)
counts_fast = create_counts_fast(id_sort)
res = euler_full()
# counts_np = np.zeros((T, N), dtype=int)
# counts_np = create_counts_np(id_sort, counts_np)

# %%

# res = eulerian_from_count(data_sort, counts_fast, N)
# M = ak.max(ak.num(res, axis=-1))
# res_padded = ak.pad_none(res, M, axis=-1)
# ma = res_padded.to_numpy()
# ma_float = ma.astype(float)


# a = np.ma.getdata(ma_float)
# a[np.ma.getmask(ma_float)] = np.nan
# a

# %%


def get_shape(a):
    np_size = [int(ak.num(a, axis=0))]
    for i in range(1, a.ndim):
        np_size.append(ak.max(ak.num(res, axis=i)))
    return np_size


def calc_memory_array(elements):
    memory_usage = elements * 8  # 8 bytes per float

    if memory_usage < 1024 * 1024:
        memory_usage_mb = memory_usage / (1024 * 1024)
        return f"Memory usage: {memory_usage_mb:.2f} Mb"
    elif memory_usage < 1024 * 1024 * 1024:
        memory_usage_gb = memory_usage / (1024 * 1024 * 1024)
        return f"Memory usage: {memory_usage_gb:.2f} Gb"
    else:
        memory_usage_tb = memory_usage / (1024 * 1024 * 1024 * 1024)
        return f"Memory usage: {memory_usage_tb:.2f} Tb"


print(calc_memory_array(np.prod(get_shape(res))))
print(calc_memory_array(ak.count(res)))
# %%
