"""
Copyright (c) 2024 MPI-M, Clara Bayley


----- CLEO -----
File: yac_data_reader.py
Project: yac_3d
Created Date: Tuesday 7th May 2024
Author: Wilton Loch (WL)
Additional Contributors:
-----
Last Modified: Tuesday 7th May 2024
Modified By: CB
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
"""

from yac import YAC, def_calendar, Calendar, Reg2dGrid, Location, Field, TimeUnit
import numpy as np

# --- Subroutines for binary data reading ---


class variable_metadata:
    def __init__(self, b0, bsize, nvar, vtype, units, scale_factor):
        self.b0 = b0
        self.bsize = bsize
        self.nvar = nvar
        self.vtype = vtype
        self.units = units
        self.scale_factor = scale_factor


def thermodynamicvar_from_binary(filename):
    file = open(filename, "rb")
    print("Opened file", filename)

    variable_metadata = read_metadata(file)

    return vector_from_binary(file, variable_metadata[0])


def vector_from_binary(file, variable_metadata):
    file.seek(variable_metadata.b0)
    data = np.array(
        np.frombuffer(file.read(8 * variable_metadata.nvar), dtype=np.float64),
        copy=True,
    )
    return data


def read_variable_metadata(file, offset):
    b0 = int.from_bytes(file.read(4), "little")
    bsize = int.from_bytes(file.read(4), "little")
    nvar = int.from_bytes(file.read(4), "little")

    vtype = file.read(1).decode("UTF-8")
    units = file.read(1).decode("UTF-8")

    scale_factor = np.frombuffer(file.read(8), dtype=np.float64)[0]

    return variable_metadata(b0, bsize, nvar, vtype, units, scale_factor)


def read_metadata(file):
    file.read(4)
    charbytes = int.from_bytes(file.read(4), "little")
    nvars = int.from_bytes(file.read(4), "little")
    mbytes_pervar = int.from_bytes(file.read(4), "little")
    file.read(charbytes)

    variable_metadata = []
    offset = 4 + charbytes

    for var_index in range(nvars):
        variable_metadata.append(read_variable_metadata(file, offset))
        offset = offset + mbytes_pervar

    return variable_metadata


# --- Start of YAC definitions ---
yac = YAC()

# --- Component definition ---
component_name = "yac_reader"
component = yac.def_comp(component_name)

def_calendar(Calendar.PROLEPTIC_GREGORIAN)

# --- Grid definition ---
lon = np.linspace(0, 2 * np.pi, 27)[:-1]
lat = np.linspace(-0.5 * np.pi, 0.5 * np.pi, 33)[1:-1]
grid = Reg2dGrid("yac_reader_grid", lon, lat)

# --- Point definitions ---
cell_centers_lon = (lon + np.pi / 27)[:-1]
cell_centers_lat = (lat + np.pi / 66)[:-1]
edge_centers_lat = []
edge_centers_lon = []

for lat_index in range(0, len(lat) * 2 - 1):
    if lat_index % 2 == 0:
        edge_centers_lon.extend(cell_centers_lon)
        edge_centers_lat.extend([lat[lat_index // 2]] * len(cell_centers_lon))
    else:
        edge_centers_lon.extend(lon)
        edge_centers_lat.extend([cell_centers_lat[(lat_index - 1) // 2]] * len(lon))

cell_centers = grid.def_points(Location.CELL, cell_centers_lon, cell_centers_lat)
edge_centers = grid.def_points_unstruct(
    Location.EDGE, edge_centers_lon, edge_centers_lat
)

# --- Field definitions ---
press = Field.create(
    "pressure", component, cell_centers, 1, "PT1M", TimeUnit.ISO_FORMAT
)
temp = Field.create(
    "temperature", component, cell_centers, 1, "PT1M", TimeUnit.ISO_FORMAT
)
qvap = Field.create("qvap", component, cell_centers, 1, "PT1M", TimeUnit.ISO_FORMAT)
qcond = Field.create("qcond", component, cell_centers, 1, "PT1M", TimeUnit.ISO_FORMAT)
vvel = Field.create("vvel", component, cell_centers, 1, "PT1M", TimeUnit.ISO_FORMAT)
hor_wind_velocities = Field.create(
    "hor_wind_velocities", component, edge_centers, 1, "PT1M", TimeUnit.ISO_FORMAT
)

# --- End of YAC definitions ---
yac.enddef()
np.set_printoptions(threshold=np.inf)

# Read binary data from files
press_values = thermodynamicvar_from_binary(
    "../build/share/yac1_dimlessthermo_press.dat"
)
temp_values = thermodynamicvar_from_binary("../build/share/yac1_dimlessthermo_temp.dat")
qvap_values = thermodynamicvar_from_binary("../build/share/yac1_dimlessthermo_qvap.dat")
qcond_values = thermodynamicvar_from_binary(
    "../build/share/yac1_dimlessthermo_qcond.dat"
)
uvel = thermodynamicvar_from_binary("../build/share/yac1_dimlessthermo_uvel.dat")
wvel = thermodynamicvar_from_binary("../build/share/yac1_dimlessthermo_wvel.dat")
vvel_values = thermodynamicvar_from_binary("../build/share/yac1_dimlessthermo_vvel.dat")

# Pack all horizontal wind velocities edge data into one united field for YAC exchange
united_edge_data = []
for timestep in range(5):
    u_timestep_offset = timestep * 2325
    w_timestep_offset = timestep * 2340
    for vertical_level in range(3):
        for lat_index in range(len(lat) * 2 - 1):
            if lat_index % 2 == 0:
                vertical_level_offset = vertical_level * 775
                lower_index = (
                    u_timestep_offset
                    + vertical_level_offset
                    + (lat_index // 2) * (len(lon) - 1)
                )
                upper_index = (
                    u_timestep_offset
                    + vertical_level_offset
                    + (lat_index // 2 + 1) * (len(lon) - 1)
                )
                united_edge_data.extend(uvel[lower_index:upper_index])
            else:
                vertical_level_offset = vertical_level * 780
                lower_index = (
                    w_timestep_offset
                    + vertical_level_offset
                    + ((lat_index - 1) // 2) * len(lon)
                )
                upper_index = (
                    w_timestep_offset
                    + vertical_level_offset
                    + ((lat_index + 1) // 2) * len(lon)
                )
                united_edge_data.extend(wvel[lower_index:upper_index])

# somewhat generic version for cell and horizontal edge coupling
for timestep in range(5):
    timestep_cell_offset = timestep * 2250
    timestep_edge_offset = timestep * 4665
    vertical_winds_lower_index = timestep * 3000
    vertical_winds_upper_index = vertical_winds_lower_index + 750

    vvel.put(vvel_values[vertical_winds_lower_index:vertical_winds_upper_index])

    for vertical_level in range(3):
        vertical_level_cell_offset = vertical_level * 750
        vertical_level_edge_offset = vertical_level * 1555

        cell_lower_index = timestep_cell_offset + vertical_level_cell_offset
        cell_upper_index = cell_lower_index + 750

        edge_lower_index = timestep_edge_offset + vertical_level_edge_offset
        edge_upper_index = edge_lower_index + 1555

        vertical_winds_lower_index = vertical_winds_upper_index
        vertical_winds_upper_index = vertical_winds_lower_index + 750

        vvel.put(vvel_values[vertical_winds_lower_index:vertical_winds_upper_index])
        press.put(press_values[cell_lower_index:cell_upper_index])
        temp.put(temp_values[cell_lower_index:cell_upper_index])
        qvap.put(qvap_values[cell_lower_index:cell_upper_index])
        qcond.put(qcond_values[cell_lower_index:cell_upper_index])
        hor_wind_velocities.put(
            np.asarray(united_edge_data[edge_lower_index:edge_upper_index])
        )

del yac
