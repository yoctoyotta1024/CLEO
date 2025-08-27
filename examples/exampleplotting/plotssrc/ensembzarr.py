"""
Copyright (c) 2025 MPI-M, Clara Bayley


----- CLEO -----
File: ensembzarr.py
Project: plotssrc
Created Date: Friday 17th November 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
File Description:
functions to write a new zarr dataset and setuptxt file for an ensemble of datasets
"""


import shutil
import sys
import numpy as np
import zarr

from ....cleopy.sdmout_src import pyzarr
from ....cleopy import editconfigfile


def write_ensemble_info(ensembsetupfile, setupfile, datasets):
    """add extra information in ensembsetupfile about
    which datasets and setuptxt file were used to make ensemble"""

    header = (
        "// ----------------------------- //"
        + "\n// --------- ensembsetupfile --------- //"
        + "\n// ----------------------------- //\n"
    )
    footer = "// ----------------------------- //"

    datasets_str = "\ndatasets in ensemble: \n     " + "\n     ".join(datasets) + "\n"

    setup_str = "\nsetup copied from: \n     " + setupfile + "\n"

    with open(ensembsetupfile, "a") as file:
        file.write(header)
        file.write(setup_str)
        file.write(datasets_str)
        file.write(footer)


def write_ensemb_setupfile(ensembsetupfile, setupfile, datasets):
    """copy setupfile to ensembsetupfile and then edit / add
    information relevant to datasets used to make ensemble"""

    shutil.copy(setupfile, ensembsetupfile)
    params = {
        "initsupers_filename": "[ensemble, see below]",
        "setup_filename": "[ensemble, see below]",
        "zarrbasedir": "[ensemble, see below]",
    }
    editconfigfile.edit_config_params(ensembsetupfile, params)

    write_ensemble_info(ensembsetupfile, setupfile, datasets)


def check_dataset_for_ensemb(dataset, refset):
    """returns xarray dataset after checking that
    time is consistent with the refset"""

    reftime = pyzarr.get_rawdataset(refset)["time"].values
    ds = pyzarr.get_rawdataset(dataset)

    try:
        time = ds["time"].values
        if np.any(time != reftime):
            print("refset: " + refset + ", dataset: " + dataset)
            raise ValueError("data for time in dataset must be same as reference")
    except KeyError:
        raise KeyError("no time data in dataset " + dataset)


def write_time_to_ensembzarr(ensembdataset, dataset):
    """create or replace time group in ensembdataset
    zarr storage with time group copied from dataset"""

    src = zarr.open(dataset)["time"]
    dest = zarr.open(ensembdataset)
    zarr.copy(src, dest, log=sys.stdout, if_exists="replace")


def ensemble_data(func, datasets, var):
    """call func on ensemble of var from datasets.
    Note: if var has dimensions [a,b, ...] then data4ensemb
    passed into func has dimensions [n, a, b, ...]
    when n is the number of datasets in the ensemble"""

    data4ensemb = []
    for dataset in datasets:
        v = pyzarr.get_rawdataset(dataset)[var].values
        data4ensemb.append(v)
    data4ensemb = np.asarray(data4ensemb)

    return func(data4ensemb)


def write_matchingarray_to_storage(arrayname, arr, refz, zattrs):
    """write array 'arr' to zarr storage under
    'arrayname' using same metadata as refz and
    copying zattrs to .zattrs file"""

    if arr.shape != refz.shape:
        err = "arr and reference for metadata must have same shape"
        raise ValueError(err)

    z = zarr.open(
        arrayname,
        mode="w",
        shape=refz.shape,
        chunks=refz.chunks,
        compressor=refz.compressor,
        dtype=refz.dtype,
        fill_value=refz.fill_value,
        filters=refz.filters,
        order=refz.order,
    )
    z[:] = arr

    shutil.copy(zattrs, arrayname + "/.zattrs")


def write_meanvars_to_ensembzarr(ensembdataset, vars4ensemb, datasets, refset):
    """write mean over ensemble of datasets for
    each var in vars4ensemb into zarr array also
    called 'var' in ensembdataset"""

    for var in vars4ensemb:
        meanname = ensembdataset + "/" + var
        meanvar = ensemble_data(lambda x: np.mean(x, axis=0), datasets, var)

        zattrs = refset + "/" + var + "/.zattrs"
        write_matchingarray_to_storage(
            meanname, meanvar, zarr.open(refset)[var], zattrs
        )


def write_stdvars_to_ensembzarr(ensembdataset, vars4ensemb, datasets, refset):
    """write standard deviation over ensemble of datasets
    for each var in vars4ensemb into zarr array also
    called 'var' in ensembdataset"""

    for var in vars4ensemb:
        stdname = ensembdataset + "/" + var + "_std"
        stdvar = ensemble_data(lambda x: np.std(x, axis=0), datasets, var)
        stdvar = stdvar / np.sqrt(len(datasets))

        zattrs = refset + "/" + var + "/.zattrs"
        write_matchingarray_to_storage(stdname, stdvar, zarr.open(refset)[var], zattrs)


def write_ensemb_zarr(ensembdataset, vars4ensemb, datasets):
    """create zarr storage called 'ensembdataset'
    for an ensemble of zarr datasets 'datasets'
    containing their time data alongside their
    mean and standard deviation for the variables listed in
    'vars4ensemb'"""

    refset = datasets[0]  # reference dataset
    for dataset in datasets:
        check_dataset_for_ensemb(dataset, refset)

    write_time_to_ensembzarr(ensembdataset, refset)

    write_meanvars_to_ensembzarr(ensembdataset, vars4ensemb, datasets, refset)

    write_stdvars_to_ensembzarr(ensembdataset, vars4ensemb, datasets, refset)


def write_ensemble(ensembdataset, ensembsetupfile, vars4ensemb, setupfile, datasets):
    write_ensemb_setupfile(ensembsetupfile, setupfile, datasets)
    write_ensemb_zarr(ensembdataset, vars4ensemb, datasets)
