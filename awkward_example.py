# %%
import awkward as ak
import numpy as np
import warnings
from tqdm import tqdm


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

id = ak.Array(
    [
        [0, 1, 2],
        [],
        [2, 3],
        [0],
    ]
)

akward_array_to_lagrange_array(
    data=data,
    dim1=time_index,
    dim2=id,
    dim1_as_index=True,
    check_indices_uniqueness=False,
)


# %%

time_unique = np.unique(time_index)
gridbox_unique = np.unique(ak.flatten(gridbox))
T = int(ak.num(time_unique, axis=0))
G = int(ak.num(gridbox_unique, axis=0))

# %%
gridbox_argsort = ak.argsort(gridbox)
gridbox_sort = gridbox[gridbox_argsort]
data_sort = data[gridbox_argsort]
id_sort = id[gridbox_argsort]

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
# res = euler2(data_sort, gridbox_sort, T, G)

from pySD.sdmout_src import pysetuptxt, pyzarr

setupfile = "/home/m/m301096/CLEO/data/output/raw/no_aerosols_collision_many_5012/clusters_301/eurec4a1d_setup.txt"
dataset = "/home/m/m301096/CLEO/data/output/raw/no_aerosols_collision_many_5012/clusters_301/eurec4a1d_sol.zarr"
# read in constants and intial setup from setup .txt file
config = pysetuptxt.get_config(setupfile, nattrs=3, isprint=False)
consts = pysetuptxt.get_consts(setupfile, isprint=False)
# Create a first simple dataset to have the coordinates for later netcdf creation
sddata = pyzarr.get_supers(dataset, consts)
# %%
# %%

# %%
time_unique = np.unique(sddata.time)
gridbox_unique = np.unique(ak.flatten(sddata.sdgbxindex))
T = int(ak.num(time_unique, axis=0))
G = int(ak.num(gridbox_unique, axis=0))

# %%
gridbox_argsort = ak.argsort(sddata.sdgbxindex)
gridbox_sort = sddata.sdgbxindex[gridbox_argsort]
data_sort = sddata.sdId[gridbox_argsort]
id_sort = sddata.sdgbxindex[gridbox_argsort]

# %%
result = euler2(data_sort, gridbox_sort, T, G)


# %%

data = ak.Array(
    [
        [10, 20, 30],
        [],
        [40, 50],
        [90],
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

id = ak.Array(
    [
        [0, 1, 2],
        [],
        [2, 3],
        [0],
    ]
)
