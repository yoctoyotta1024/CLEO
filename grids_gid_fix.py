"""
Copyright (c) 2026 MPI-M, Clara Bayley


----- CLEO -----
File: grids_gid_fix.py
Project: CLEO
Created Date: Thursday 1st January 1970
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
"""


import netCDF4
import shutil
import numpy as np

src = "/work/mh0731/m300950/icon/icon/experiments/bubble_cleo/bin/all_grids_debug.nc"  # expand ~ yourself
dst = src.replace(".nc", "_gid.nc")
shutil.copy(src, dst)

name = "cleo_cartesian_grid"
with netCDF4.Dataset(dst, "a") as ds:
    n = len(ds.dimensions[f"nc_{name}"])
    v = ds.createVariable(f"{name}.gid", "i4", (f"nc_{name}",))
    v[:] = np.arange(n)  # cell index as the "global id"
