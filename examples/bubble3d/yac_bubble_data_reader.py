"""
Copyright (c) 2024 MPI-M, Clara Bayley


----- CLEO -----
File: yac_bubble_data_reader.py
Project: bubble3d
Created Date: Friday 19th July 2024
Author: Wilton Loch (WL)
Additional Contributors: Clara Bayley (CB)
-----
Last Modified: Friday 16th August 2024
Modified By: CB
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
Python file for yac for a one-way coupling between ICON's
output data and CLEO (e.g. for the bubble test case)
"""


#!/usr/bin/env python3

from yac import YAC, UnstructuredGrid, Field, Location, Calendar, TimeUnit, def_calendar
from netCDF4 import Dataset
import numpy as np
import sys


def convert_seconds_to_isodate(seconds):
    """converts seconds interval into an ISO 8601 format string"""
    import isodate
    from datetime import timedelta

    duration = timedelta(seconds=seconds)

    return isodate.duration_isoformat(duration)


def add_seconds_to_isodate(iso, seconds):
    """adds seconds interval to ISO 8601 isodate given as a string
    and returns a string for the new ISO 8601 isodate."""

    from datetime import datetime, timedelta

    dt = datetime.fromisoformat(iso.replace("Z", "+00:00"))
    new_dt = dt + timedelta(seconds=seconds)
    new_iso = new_dt.isoformat().replace("+00:00", "Z")

    return new_iso


def map_vertices(array):
    vertex_dict = {}
    vertex_mapping = 0
    vertex_mapping_array = []

    for element in array:
        element_tuple = tuple(
            element
        )  # Convert the array to a tuple to make it hashable
        if element_tuple not in vertex_dict:
            vertex_dict[element_tuple] = vertex_mapping
            vertex_mapping += 1
        vertex_mapping_array.append(vertex_dict[element_tuple])

    return vertex_mapping_array


def create_yac_unstructured_grid(grid_filename, grid_name):
    """create unstructured grid for YAC from ICON netcdf file.
    Note different .nc files may have different names for variables,
    e.g. "nv" <=> "vertices"
    "cell" <=> "ncells"
    "clon_vertices" <=> "clon_bnds"
    "clat_vertices" <=> "clat_bnds"
    """
    # Open the NetCDF file
    dataset = Dataset(grid_filename, "r")

    # Read the variables
    nv = dataset.dimensions["nv"].size
    no_cells = dataset.dimensions["cell"].size

    vertices = []
    for cell in range(len(dataset["clon_vertices"])):
        for vertex in range(3):
            vertices.append(
                [
                    dataset["clon_vertices"][cell][vertex],
                    dataset["clat_vertices"][cell][vertex],
                ]
            )

    cell_vertex_indices = map_vertices(vertices)

    grid = UnstructuredGrid(
        grid_name,
        np.ones(no_cells) * nv,
        dataset["clon_vertices"][:, :].flatten(),
        dataset["clat_vertices"][:, :].flatten(),
        cell_vertex_indices,
    )

    cell_points = grid.def_points(Location.CELL, dataset["clon"][:], dataset["clat"][:])
    grid.cell_points = cell_points

    return grid


def prepare_data_for_yac(source):
    vertical_levels = len(source[:, 1])
    horizontal_size = len(source[1, :])
    target = np.empty(shape=(vertical_levels, 1, horizontal_size))
    for level in range(vertical_levels):
        target[level, 0, :] = source[level, :]
    return target


grid_filename = sys.argv[1]
data_filename = sys.argv[2]
grid_name = sys.argv[3]  # must match CLEO (see yac_cartesian_dynamics.cpp)
DATATSTEP = float(sys.argv[4])  # must match ICON data file [seconds]
COUPLTSTEP = float(sys.argv[5])  # must match CLEO config file [seconds]
T_END = float(sys.argv[6])  # must match CLEO config file [seconds]
num_vertical_levels = int(sys.argv[7])  # must match CLEO gridfile

msg = (
    "--- INPUT ARGS ---"
    + "\ngrid_filename: "
    + grid_filename
    + "\ndata_filename: "
    + data_filename
    + "\ngrid_name: "
    + grid_name
    + "\nDATASTEP: "
    + str(DATATSTEP)
    + "\nCOUPLTSTEP: "
    + str(COUPLTSTEP)
    + "\nT_END: "
    + str(T_END)
    + "\nnum_vertical_levels: "
    + str(num_vertical_levels)
    + "\n--- --- --- --- --- ---"
)
print(msg)

assert (
    DATATSTEP <= COUPLTSTEP
), "COUPLTSTEP [s] must be greater than or equal to DATATSTEP [s]"
assert (
    COUPLTSTEP % DATATSTEP == 0.0
), "COUPLTSTEP [s] must be integer multiple of DATATSTEP [s]"

yac = YAC()

def_calendar(Calendar.PROLEPTIC_GREGORIAN)
iso_start = "2008-08-01T00:00:00Z"
iso_end = add_seconds_to_isodate(iso_start, T_END)
yac.def_datetime(iso_start, iso_end)

component_name = "atm"
component = yac.def_comp(component_name)
grid = create_yac_unstructured_grid(grid_filename)

# --- Field definitions ---
coupling_tstep = convert_seconds_to_isodate(COUPLTSTEP)

press = Field.create(
    "pressure",
    component,
    grid.cell_points,
    num_vertical_levels,
    coupling_tstep,
    TimeUnit.ISO_FORMAT,
)
temp = Field.create(
    "temperature",
    component,
    grid.cell_points,
    num_vertical_levels,
    coupling_tstep,
    TimeUnit.ISO_FORMAT,
)
qvap = Field.create(
    "qvap",
    component,
    grid.cell_points,
    num_vertical_levels,
    coupling_tstep,
    TimeUnit.ISO_FORMAT,
)
qcond = Field.create(
    "qcond",
    component,
    grid.cell_points,
    num_vertical_levels,
    coupling_tstep,
    TimeUnit.ISO_FORMAT,
)
eastward_wind = Field.create(
    "eastward_wind",
    component,
    grid.cell_points,
    num_vertical_levels,
    coupling_tstep,
    TimeUnit.ISO_FORMAT,
)
northward_wind = Field.create(
    "northward_wind",
    component,
    grid.cell_points,
    num_vertical_levels,
    coupling_tstep,
    TimeUnit.ISO_FORMAT,
)
vertical_wind = Field.create(
    "vertical_wind",
    component,
    grid.cell_points,
    num_vertical_levels + 1,
    coupling_tstep,
    TimeUnit.ISO_FORMAT,
)

yac.enddef()

dataset = Dataset(data_filename)

datasteps_per_coupling_step = COUPLTSTEP / DATATSTEP
num_couplingsteps = int(np.floor(T_END / COUPLTSTEP) + 1)

for coupling_step in range(num_couplingsteps):
    timeindex = (
        coupling_step * datasteps_per_coupling_step
    )  # index along time axis of data to "put"
    temp.put(prepare_data_for_yac(dataset["ta"][timeindex, 0:num_vertical_levels, :]))
    press.put(
        prepare_data_for_yac(dataset["pfull"][timeindex, 0:num_vertical_levels, :])
    )
    qvap.put(prepare_data_for_yac(dataset["hus"][timeindex, 0:num_vertical_levels, :]))
    qcond.put(prepare_data_for_yac(dataset["clw"][timeindex, 0:num_vertical_levels, :]))
    vertical_wind.put(
        prepare_data_for_yac(dataset["wa"][timeindex, 0 : num_vertical_levels + 1, :])
    )
    eastward_wind.put(
        prepare_data_for_yac(dataset["ua"][timeindex, 0:num_vertical_levels, :])
    )
    northward_wind.put(
        prepare_data_for_yac(dataset["va"][timeindex, 0:num_vertical_levels, :])
    )
