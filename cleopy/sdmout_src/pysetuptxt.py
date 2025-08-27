"""
Copyright (c) 2024 MPI-M, Clara Bayley


----- CLEO -----
File: pysetuptxt.py
Project: sdmout_src
Created Date: Tuesday 24th October 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
functions for reading setup.txt file
output alongside zarr storage
"""

from .. import cxx2py


def get_consts(setuptxt, isprint=True):
    """returns dictionary of constants
    read from from setup.txt file"""

    return consts_dict(setuptxt, isprint=isprint)


def get_config(setuptxt, nattrs=3, isprint=True):
    """returns dictionary of configuration parameters
    read from from setup.txt file"""

    return config_dict(setuptxt, nattrs=nattrs, isprint=isprint)


def consts_dict(setuptxt, isprint):
    consts = cxx2py.read_cxxconsts_into_floats(setuptxt)
    consts.update(cxx2py.derive_more_floats(consts))

    if isprint:
        cxx2py.print_dict_statement(setuptxt, "consts", consts)

    return consts


def config_dict(setuptxt, nattrs, isprint):
    config = read_configparams_fromsetuptxt_into_floats(setuptxt)
    config["numSDattrs"] = config["nspacedims"] + nattrs
    config["ntime"] = round(config["T_END"] / config["OBSTSTEP"]) + 1

    if isprint:
        cxx2py.print_dict_statement(setuptxt, "config", config)

    return config


def read_configparams_fromsetuptxt_into_floats(filename):
    """returns dictionary of value: float from
    values assigned in a config .txt file.
    Also returns dictionary of notfloats
    for values that couldn't be converted."""

    floats = {}
    notfloats = {}
    with open(filename) as file:
        rlines = []
        filelines = file.readlines()
        for line in filelines:
            if (line[0] != "#") and (line[0] != "/") and (":" in line):
                goodline = cxx2py.remove_excess_line(line)
                rlines.append(goodline)

        for line in rlines:
            ind = line.find(":")
            name = line[:ind]
            value = line[ind + 1 :]

            try:
                floats[name] = float(value)
            except ValueError:
                notfloats[name] = value

    try:
        floats["nspacedims"] = int(floats["nspacedims"])  # no spatial coords to SDs
    except KeyError as e:
        print("Warning, ignoring error", e)
        pass

    return floats
