'''
----- CLEO -----
File: pysetuptxt.py
Project: sdmout_src
Created Date: Tuesday 24th October 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
Last Modified: Wednesday 17th April 2024
Modified By: CB
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
Copyright (c) 2023 MPI-M, Clara Bayley
-----
File Description:
functions for reading setup.txt file
output alongside zarr storage
'''

from .. import cxx2py, readconfigfile

def get_config(setuptxt, nattrs=3, isprint=True):
  '''returns dictionary of configuration parameters
  read from from setup.txt file '''

  return config_dict(setuptxt, nattrs=nattrs, isprint=isprint)

def get_consts(setuptxt, isprint=True):
  '''returns dictionary of constants
  read from from setup.txt file '''

  return consts_dict(setuptxt, isprint=isprint)

def config_dict(setuptxt, nattrs, isprint):

  config = readconfigfile.read_configparams_into_floats(setuptxt)
  config["numSDattrs"] = config["nspacedims"] + nattrs
  config["ntime"] = round(config["T_END"]/config["OBSTSTEP"])+1

  if isprint:
    cxx2py.print_dict_statement(setuptxt, "config", config)

  return config

def consts_dict(setuptxt, isprint):
  consts = cxx2py.read_cxxconsts_into_floats(setuptxt)
  consts.update(cxx2py.derive_more_floats(consts))

  if isprint:
    cxx2py.print_dict_statement(setuptxt, "consts", consts)

  return consts
