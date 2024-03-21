from yac import *
import numpy as np

class variable_metadata:
    def __init__(self, b0, bsize, nvar, vtype, units, scale_factor):
        self.b0           = b0
        self.bsize        = bsize
        self.nvar         = nvar
        self.vtype        = vtype
        self.units        = units
        self.scale_factor = scale_factor

def thermodynamicvar_from_binary(filename):
    file = open(filename, 'rb')
    print("Opened file", filename)

    variable_metadata = read_metadata(file)

    return vector_from_binary(file, variable_metadata[0])

def vector_from_binary(file, variable_metadata):
    file.seek(variable_metadata.b0)
    data = np.array(np.frombuffer(file.read(8 * variable_metadata.nvar), dtype = np.float64), copy=True)
    return data

def read_variable_metadata(file, offset):
    b0    = int.from_bytes(file.read(4), "little")
    bsize = int.from_bytes(file.read(4), "little")
    nvar  = int.from_bytes(file.read(4), "little")

    vtype = file.read(1).decode("UTF-8")
    units = file.read(1).decode("UTF-8")

    scale_factor = np.frombuffer(file.read(8), dtype = np.float64)[0]

    return variable_metadata(b0, bsize, nvar, vtype, units, scale_factor)

def read_metadata(file):
    d0byte        = int.from_bytes(file.read(4), "little")
    charbytes     = int.from_bytes(file.read(4), "little")
    nvars         = int.from_bytes(file.read(4), "little")
    mbytes_pervar = int.from_bytes(file.read(4), "little")

    global_meta_string = file.read(charbytes).decode("UTF-8")

    variable_metadata = []
    offset = 4 + charbytes

    for var_index in range(nvars):
        variable_metadata.append(read_variable_metadata(file, offset))
        offset = offset + mbytes_pervar

    return variable_metadata

yac = YAC()

component_name = "yac_reader"
component = yac.def_comp(component_name)

def_calendar(Calendar.PROLEPTIC_GREGORIAN)

lon = np.linspace(0,2*np.pi,32)[:-1]
lat = np.linspace(-0.5*np.pi,0.5*np.pi, 33)[1:-1]
cell_centers_lon = (lon + np.pi/32)[:-1]
cell_centers_lat = (lat + np.pi/66)[:-1]

grid = Reg2dGrid(f"yac_reader_grid", lon, lat)
cell_centers = grid.def_points(Location.CELL, cell_centers_lon, cell_centers_lat)

press = Field.create("pressure", component, cell_centers, 1, "PT1M", TimeUnit.ISO_FORMAT)
temp  = Field.create("temperature", component, cell_centers, 1, "PT1M", TimeUnit.ISO_FORMAT)
qvap  = Field.create("qvap", component, cell_centers, 1, "PT1M", TimeUnit.ISO_FORMAT)
qcond = Field.create("qcond", component, cell_centers, 1, "PT1M", TimeUnit.ISO_FORMAT)

yac.enddef()
np.set_printoptions(threshold=np.inf)

press_values = thermodynamicvar_from_binary("../build/share/df2d_dimlessthermo_press.dat")
temp_values = thermodynamicvar_from_binary("../build/share/df2d_dimlessthermo_temp.dat")
qvap_values = thermodynamicvar_from_binary("../build/share/df2d_dimlessthermo_qvap.dat")
qcond_values = thermodynamicvar_from_binary("../build/share/df2d_dimlessthermo_qcond.dat")

press.put(press_values[0:900])
temp.put(temp_values[0:900])
qvap.put(qvap_values[0:900])
qcond.put(qcond_values[0:900])

press.put(press_values[900:1800])
temp.put(temp_values[900:1800])
qvap.put(qvap_values[900:1800])
qcond.put(qcond_values[900:1800])

press.put(press_values[1800:2700])
temp.put(temp_values[1800:2700])
qvap.put(qvap_values[1800:2700])
qcond.put(qcond_values[1800:2700])

del yac
