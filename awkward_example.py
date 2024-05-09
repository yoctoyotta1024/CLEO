# %%
import awkward as ak
import numpy as np
from tqdm import tqdm

# %%
data = ak.Array(
    [
        [10, 20, 30],
        [],
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
        [0, 1, 0],
        [],
        [1, 2],
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
        [0, 1, 2],
        [7],
        [2, 3],
        [7],
    ]
)


# %%

time_unique = np.unique(time_index)
gridbox_unique = np.unique(ak.flatten(gridbox))
id_unique = np.unique(ak.flatten(id))
T = int(ak.num(time_unique, axis=0))
G = int(ak.num(gridbox_unique, axis=0))
N = int(ak.max(id_unique, axis=0)) + 1


gridbox_argsort = ak.argsort(gridbox)
gridbox_argsort = ak.argsort(id)
gridbox_sort = gridbox[gridbox_argsort]
data_sort = data[gridbox_argsort]
id_sort = id[gridbox_argsort]


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


# %%


def create_counts_fast(data, T, G):
    maximum = ak.max(data) + 1
    modified = data + np.arange(T) * maximum
    modified = ak.flatten(modified)

    counts_flat = np.bincount(modified)
    counts_flat = ak.unflatten(counts_flat, G)
    return counts_flat


eulerian_from_count(data_sort, create_counts_fast(id_sort, T, N), N)
# unique = ak.unflatten(unique, num)
# counts = ak.unflatten(counts, num)

# %%


def cut_unique_values(array, data):
    if len(array) != 0:
        unique, counts = np.unique(array, return_counts=True, axis=0)
        return ak.unflatten(data, counts)
    else:
        return array


# %%
cut_unique_values(gridbox_sort[0], data_sort[0])


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

from pySD.sdmout_src import pysetuptxt, pyzarr
import numpy as np

setupfile = "/home/m/m301096/CLEO/data/output/raw/no_aerosols_collision_many_5012/clusters_301/eurec4a1d_setup.txt"
dataset = "/home/m/m301096/CLEO/data/output/raw/no_aerosols_collision_many_5012/clusters_301/eurec4a1d_sol.zarr"
# read in constants and intial setup from setup .txt file
config = pysetuptxt.get_config(setupfile, nattrs=3, isprint=False)
consts = pysetuptxt.get_consts(setupfile, isprint=False)
# Create a first simple dataset to have the coordinates for later netcdf creation
sddata = pyzarr.get_supers(dataset, consts)

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
counts_fast = create_counts_fast(id_sort, T, N)
counts_np = np.zeros((T, N), dtype=int)
counts_np = create_counts_np(id_sort, counts_np)

# %%

res = eulerian_from_count(data_sort, counts_fast, N)
M = ak.max(ak.num(res, axis=-1))
res_padded = ak.pad_none(res, M, axis=-1)
ma = res_padded.to_numpy()
ma_float = ma.astype(float)


a = np.ma.getdata(ma_float)
a[np.ma.getmask(ma_float)] = np.nan
a
