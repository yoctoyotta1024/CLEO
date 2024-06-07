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

    grid_name = "bubble_grid"
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


data_filename = "aes_bubble_atm_3d_ml_20080801T000000Z.nc"
grid_filename = "aes_bubble_atm_cgrid_ml.nc"

yac = YAC()

component_name = "icon_data_reader"
component = yac.def_comp(component_name)

def_calendar(Calendar.PROLEPTIC_GREGORIAN)
grid = create_yac_unstructured_grid(grid_filename)

# --- Field definitions ---
press = Field.create(
    "pressure", component, grid.cell_points, 1, "PT1M", TimeUnit.ISO_FORMAT
)
temp = Field.create(
    "temperature", component, grid.cell_points, 1, "PT1M", TimeUnit.ISO_FORMAT
)
qvap = Field.create("qvap", component, grid.cell_points, 1, "PT1M", TimeUnit.ISO_FORMAT)
qcond = Field.create(
    "qcond", component, grid.cell_points, 1, "PT1M", TimeUnit.ISO_FORMAT
)
eastward_wind = Field.create(
    "eastward_wind", component, grid.cell_points, 1, "PT1M", TimeUnit.ISO_FORMAT
)
northward_wind = Field.create(
    "northward_wind", component, grid.cell_points, 1, "PT1M", TimeUnit.ISO_FORMAT
)
vvel = Field.create("vvel", component, grid.cell_points, 1, "PT1M", TimeUnit.ISO_FORMAT)

yac.enddef()

dataset = Dataset(data_filename)

for step in range(5):
    print("ICON data reader started sending step", step)
    vvel.put(dataset["wa"][step, 0, :])
    for height in range(3):
        temp.put(dataset["ta"][step, height, :])
        press.put(dataset["pfull"][step, height, :])
        qvap.put(dataset["hus"][step, height, :])
        qcond.put(dataset["clw"][step, height, :])
        eastward_wind.put(dataset["ua"][step, height, :])
        northward_wind.put(dataset["va"][step, height, :])
        vvel.put(dataset["wa"][step, height + 1, :])
    print("ICON data reader finished sending step", step)
