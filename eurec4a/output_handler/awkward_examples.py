# %%
import numpy as np
import awkward as ak
import sys
import numpy as np
from pathlib import Path

# from plotssrc import pltsds, pltmoms, animations
from pySD.sdmout_src import *
from pySD.initsuperdropsbinary_src import *
from pySD.sdmout_src.supersdata import SupersData 

from sdm_eurec4a.visulization import set_custom_rcParams
set_custom_rcParams()

path2CLEO = Path("/home/m/m301096/CLEO")
rawdirectory = path2CLEO / "data/output/raw/rain"

sys.path.append(path2CLEO)  # for imports from pySD package
# sys.path.append(path2CLEO / "examples/exampleplotting") # for imports from example plotting package

# use paths to files
path2build = path2CLEO / "build"
configfile = path2CLEO / "eurec4a/experiment_03/src/config/rain1d_config.txt"
# yaml_config_file = path2sdm_eurec4a / "data/model/input/example_input_18.yaml"

constsfile    = path2CLEO / "libs/cleoconstants.hpp"
# path and file names for plotting results
setupfile     = path2CLEO / "data/output/raw/rain/clusters_18/rain1d_setup.txt"
dataset       = path2CLEO / "data/output/raw/rain/clusters_18/rain1d_sol.zarr"

# read in constants and intial setup from setup .txt file
config = pysetuptxt.get_config(setupfile, nattrs=3, isprint=False)
consts = pysetuptxt.get_consts(setupfile, isprint=False)

sddata = SupersData(
    dataset = str(dataset),
    consts = consts
)
lagrange = sddata.to_Dataset()

# %%
def ak_apply(array, np_func, kwargs) :
    counts = ak.num(array)

    flat = ak.flatten(array)
    flat.type.show()
    flat = np_func(flat, **kwargs)

    return ak.unflatten(array=flat, counts=counts) 

# %%


res = akward_array_to_lagrange_array(sddata["xi"], sddata.time, sddata["sdgbxindex"], check_indices_uniqueness=True)

# %%
# create the output dimensions of the numpy array which are necessary to store the data


data = ak.Array([
        [10, 20, 30],
        [],
        [40, 50],
    ])
time = ak.Array([
        10.0,
        20.0,
        45.0,
    ])
id = ak.Array([
        [0, 1, 2],
        [],
        [2, 1],
    ])
x0 = ak.Array([
        [0, 0, 1],
        [],
        [1, 1],
    ])
x1 = ak.Array([
        [1, 0, 0],
        [],
        [0, 0],
    ])
x2 = ak.Array([
        [0, 0, 0],
        [],
        [0, 0],
    ])

# data = sddata["xi"]
# time = sddata.time
# id = sddata["sdId"]
# x0 = sddata["coord1"]
# x1 = sddata["coord2"]
# x2 = sddata["coord3"]
# x2 = ak_apply(x2, np.digitize, {"bins" : np.arange(0, 1200, 20)})

reduction_dim = id
dim1 = time
dim1_as_index = False
if dim1_as_index is False:
    T = int(ak.num(dim1, axis = 0))
    leading_dim = np.arange(T)
else :
    leading_dim = dim1
    T = int(ak.max(dim1) + 1)

# The superdroplets are identified by their id.
# Use the maximum value of the superdroplet index
# The ids start with id "0", so the maximum id is the number of superdroplets - 1!
S = int(ak.max(reduction_dim) + 1)

def get_index(dim : ak.Array, leading_dim = leading_dim, reduction_dim = reduction_dim) :
    if ak.count(dim) > 0:
        L = int(ak.max(dim) + 1)
        t, l = ak.unzip(ak.flatten(ak.cartesian((leading_dim, dim), axis = 1)))
    else :
        L = 1
        t, l = ak.unzip(ak.flatten(ak.cartesian((leading_dim, reduction_dim), axis = 1)))
        l = l * 0
    return L, t, l

M, t0, i = get_index(x0)
N, t1, j = get_index(x1)
O, t2, k = get_index(x2)

# check if the resulting array is sparse and inform the User

# The list of tuples is then unzipped, to seperate the time and superdroplet indeices into two arrays
t_red, s = ak.unzip(ak.flatten(ak.cartesian((leading_dim, reduction_dim), axis = 1)))
if np.all((t0 == t1) & (t0 == t2)) == True:
    t = t0
else:
    raise ValueError("The time indices are not the same for all dimensions.")

# check if the resulting array is sparse and inform the User
filled_percentage = ak.count(reduction_dim) / (T * M * N * O * S) * 100
print(f"{filled_percentage:.2f} % of the regular array is filled with values.")
    
result_numpy = np.full((T, M, N, O, S), np.nan)
result_numpy[t, i, j, k, s] = ak.flatten(data)

result = np.sum(result_numpy, axis = -1, where=~np.isnan(result_numpy))
result[(~np.isnan(result_numpy)).sum(axis = -1) > 0] = np.nan
100*(1 - (np.sum(np.isnan(result)) / result.size))
# da = xr.DataArray(
#     result_numpy,
#     dims = ["time", "x0", "x1", "x2", "id"],
#     coords = {
#         "time" : time.to_list(),
#         "x0" : np.arange(M),
#         "x1" : np.arange(N),
#         "x2" : np.arange(O),
#         "id" : np.arange(S),
#         },
# )
# da.mean("id")
# Check if the indice tuples are unique.
# For this, one can simply compute the index value they represented in a raveled array.
# so i * N + j should be unique for all i, j
# flattened_index = i * N + j
# if check_indices_uniqueness is True :
#     if not ak.count(i * N + j) == ak.count(np.unique(i * N + j)) :
#         raise ValueError(
#             "The indice tuples are not unique.\n"\
#             +"This would lead to overwriting values in the numpy array.\n"\
#             +"The reason might be, that the time indices aren't unique along axis 0 already.\n"
#             )
# else :
#     warnings.warn("The uniqueness of the indices is not checked. This might lead to overwriting values in the numpy array.")

# result_numpy = np.empty((T, N)) * np.nan
# result_numpy[i, j] = ak.flatten(data)
# return result_numpy



# result = akward_array_to_lagrange_array(sddata["xi"], sddata.time, sddata["sdId"])
# %%
# t = np.arange(len(sddata.time))
# i = sddata.sdId
# a = sddata.coord3
# counts = ak.num(a)

# T = int(ak.num(t, axis = 0))

# N = int(ak.max(i) + 1)

# a_pad = ak.pad_none(a, target = 3)
# a_pad = ak.fill_none(a_pad, np.nan)
# counts_pad = ak.num(a_pad)
# res_np = np.arange(T * N).reshape(T, N) * np.nan
# res_np
# # %%
# ids = ak.cartesian((t, i))
# # ids_com = ak.combinations((t, i), n = 2)

# ids_flat = ak.flatten(ids, axis = -1)
# # t_idx = ids_flat[..., :0]
# # sd_idx = ids_flat[..., :1]
# id_np = np.array(ids_flat.to_list()).T
# # %%
# res_np[id_np[0], id_np[1]] = ak.flatten(a)
# res_np











# #  DOES NOT WORK BECAUSE NO INFORMATION WHERE TIME AND SDID ARE RELATED
# # %%
# ds = sddata.ds
# ds = ds.drop_vars(["massmom0", "massmom1", "massmom2", "nsupers", "rgd_totnsupers"])
# ds
# # %%
# g = ds.groupby("sdId")
# # %%

# # %%

# subds = g[0]
# # def goto2D(dataset: xr.Dataset) -> xr.Dataset :
# #%%
# subds = g[22]

# dim1 = "time"
# dim2 = "sdId"
# dim2_old = "sdId"
# def swap_dims(subds: xr.Dataset, dim1: str, dim2_old: str) -> xr.Dataset:

#     dim2 = "foo"
#     subds = subds.rename({dim2_old : dim2})
#     dim2_unique = np.unique(subds[dim2].data)[0]
#     print(dim2_unique)
#     subds[dim1] = subds[dim1][:len(subds[dim2])]
#     dim_temp = subds[dim1].assign_coords({dim2 : (dim1, subds[dim2].data)})
#     dim_temp = dim_temp.set_xindex((dim2))
#     dim_temp = dim_temp.swap_dims({dim1 : dim2})
#     dim_temp = dim_temp.drop_vars(dim1)
#     subds[dim1] = dim_temp
#     subds = subds.set_xindex((dim1))
#     subds = subds.swap_dims({dim2 : dim1})
#     res = subds.drop_vars(dim2)
# # try :
# #     # res = res.assign_coords({dim2_old : (dim2_old, [dim2_unique])})
# #     pass
# # except ValueError:
# #     print(res[dim2_old])
# # res = res.expand_dims(dim2_old)
# # res[dim2_old]  = xr.DataArray([dim2_unique], dims = [dim2_old])
# # res = res.swap_dims({dim2 : dim1})

#     # return res

# # swap_dims(subds, dim1, dim2)
# # %%
# l = []
# for key, item in g:
#     l.append(swap_dims(item, dim1, dim2))
# # %%

# %%
