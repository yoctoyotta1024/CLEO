#!/usr/bin/env python3

from yac import YAC, UnstructuredGrid, Field, Location, Calendar, TimeUnit, def_calendar
from netCDF4 import Dataset
import numpy as np


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


def create_yac_unstructured_grid(grid_filname):
    # Open the NetCDF file
    dataset = Dataset(grid_filname, "r")

    # Read the variables
    nv = dataset.dimensions["vertices"].size
    no_cells = dataset.dimensions["ncells"].size

    vertices = []
    for cell in range(len(dataset["clon_bnds"])):
        for vertex in range(3):
            vertices.append(
                [dataset["clon_bnds"][cell][vertex], dataset["clat_bnds"][cell][vertex]]
            )

    cell_vertex_indices = map_vertices(vertices)

    grid_name = "icon_atmos_grid"
    grid = UnstructuredGrid(
        grid_name,
        np.ones(no_cells) * nv,
        dataset["clon_bnds"][:, :].flatten(),
        dataset["clat_bnds"][:, :].flatten(),
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


data_filename = "aes_bubble_atm_3d_ml_20080801T000000Z.nc"
grid_filename = "aes_bubble_atm_cgrid_ml.nc"

yac = YAC()

def_calendar(Calendar.PROLEPTIC_GREGORIAN)
yac.def_datetime("2008-08-01T00:00:00Z", "2008-08-01T02:00:00Z")

component_name = "atm"
component = yac.def_comp(component_name)
grid = create_yac_unstructured_grid(grid_filename)

# --- Field definitions ---
press = Field.create(
    "pressure", component, grid.cell_points, 25, "PT30M", TimeUnit.ISO_FORMAT
)
temp = Field.create(
    "temperature", component, grid.cell_points, 25, "PT30M", TimeUnit.ISO_FORMAT
)
qvap = Field.create(
    "qvap", component, grid.cell_points, 25, "PT30M", TimeUnit.ISO_FORMAT
)
qcond = Field.create(
    "qcond", component, grid.cell_points, 25, "PT30M", TimeUnit.ISO_FORMAT
)
eastward_wind = Field.create(
    "eastward_wind", component, grid.cell_points, 25, "PT30M", TimeUnit.ISO_FORMAT
)
northward_wind = Field.create(
    "northward_wind", component, grid.cell_points, 25, "PT30M", TimeUnit.ISO_FORMAT
)
vertical_wind = Field.create(
    "vertical_wind", component, grid.cell_points, 26, "PT30M", TimeUnit.ISO_FORMAT
)

yac.enddef()

dataset = Dataset(data_filename)

for step in range(5):
    temp.put(prepare_data_for_yac(dataset["ta"][step, 0:25, :]))
    press.put(prepare_data_for_yac(dataset["pfull"][step, 0:25, :]))
    qvap.put(prepare_data_for_yac(dataset["hus"][step, 0:25, :]))
    qcond.put(prepare_data_for_yac(dataset["clw"][step, 0:25, :]))
    vertical_wind.put(prepare_data_for_yac(dataset["wa"][step, 0:26, :]))
    eastward_wind.put(prepare_data_for_yac(dataset["ua"][step, 0:25, :]))
    northward_wind.put(prepare_data_for_yac(dataset["va"][step, 0:25, :]))
