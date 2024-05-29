"""
----- CLEO -----
File: rainshaft1d.py
Project: rainshaft1d
Created Date: Friday 17th November 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
Last Modified: Wednesday 17th January 2024
Modified By: CB
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
Copyright (c) 2023 MPI-M, Clara Bayley
-----
File Description:
Script compiles and runs CLEO rain1D to create the
data and plots precipitation example given constant
1-D rainshaft thermodynamics read from a file
"""

# %%
import os
import sys
import numpy as np
import random
import yaml
from pathlib import Path
from io import StringIO

print(f"Enviroment: {sys.prefix}")

path2CLEO = Path(sys.argv[1])
path2build = Path(sys.argv[2])
path2home = Path(sys.argv[3])
configfile = Path(sys.argv[4])
cloud_observation_filepath = Path(sys.argv[5])
rawdirectory = Path(sys.argv[6])

if "sdm_pysd_env312" not in sys.prefix:
    sys.path.append(str(path2CLEO))  # for imports from pySD package
    sys.path.append(
        str(path2CLEO / "examples/exampleplotting/")
    )  # for imports from example plotting package


rawdirectory.mkdir(exist_ok=True, parents=True)

# %%

with open(cloud_observation_filepath, "r") as f:
    cloud_observation_config = yaml.safe_load(f)


from pySD.sdmout_src import (
    pyzarr,
    pysetuptxt,
)
from pySD import editconfigfile
from pySD.gbxboundariesbinary_src import read_gbxboundaries as rgrid
from pySD.gbxboundariesbinary_src import create_gbxboundaries as cgrid
from pySD.initsuperdropsbinary_src import (
    crdgens,
    probdists,
    rgens,
    attrsgen,
)
from pySD.initsuperdropsbinary_src import create_initsuperdrops as csupers
from pySD.initsuperdropsbinary_src import read_initsuperdrops as rsupers
from pySD.thermobinary_src import thermogen
from pySD.thermobinary_src import create_thermodynamics as cthermo
from pySD.thermobinary_src import read_thermodynamics as rthermo


class Capturing(list):
    """
    Context manager for capturing stdout from print statements.
    https://stackoverflow.com/a/16571630/16372843
    """

    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = self._stringio = StringIO()
        return self

    def __exit__(self, *args):
        self.extend(self._stringio.getvalue().splitlines())
        del self._stringio  # free up some memory
        sys.stdout = self._stdout


# %%
### ---------------------------------------------------------------- ###
### ----------------------- INPUT PARAMETERS ----------------------- ###
### ---------------------------------------------------------------- ###
### --- essential paths and filenames --- ###
# path and filenames for creating initial SD conditions

constsfile = path2CLEO / "libs/cleoconstants.hpp"
binpath = path2build / "bin/"
sharepath = path2build / "share/"

# create the raw directory for the individual cloud observation.
# Here the output data will be stored.
identification_type = cloud_observation_config["cloud"]["identification_type"]
cloud_id = cloud_observation_config["cloud"]["cloud_id"]
rawdirectory_individual = rawdirectory / f"{identification_type}_{cloud_id}"
rawdirectory_individual.mkdir(exist_ok=True)


# CREATE A CONFIG FILE TO BE UPDATED
updated_configfile = dict()

# path and file names for output data
updated_configfile["outputdata"] = dict(
    setup_filename=str(rawdirectory_individual / "eurec4a1d_setup.txt"),
    stats_filename=str(rawdirectory_individual / "eurec4a1d_stats.txt"),
    # setup_filename = str(binpath / "eurec4a1d_setup.txt"),
    # stats_filename = str(binpath / "eurec4a1d_stats.txt"),
    zarrbasedir=str(rawdirectory_individual / "eurec4a1d_sol.zarr"),
    maxchunk=2500000,
)

setupfile = updated_configfile["outputdata"]["setup_filename"]
dataset = updated_configfile["outputdata"]["zarrbasedir"]


updated_configfile["coupled_dynamics"] = dict(
    type="fromfile",  # type of coupled dynamics to configure
    press=str(
        path2build / "share/eurec4a1d_dimlessthermo_press.dat"
    ),  # binary filename for pressure
    temp=str(
        path2build / "share/eurec4a1d_dimlessthermo_temp.dat"
    ),  # binary filename for temperature
    qvap=str(
        path2build / "share/eurec4a1d_dimlessthermo_qvap.dat"
    ),  # binary filename for vapour mixing ratio
    qcond=str(
        path2build / "share/eurec4a1d_dimlessthermo_qcond.dat"
    ),  # binary filename for liquid mixing ratio
    wvel=str(
        path2build / "share/eurec4a1d_dimlessthermo_wvel.dat"
    ),  # binary filename for vertical (coord3) velocity
    thermo=str(
        path2build / "share/eurec4a1d_dimlessthermo.dat"
    ),  # binary filename for thermodynamic profiles
)

updated_configfile["inputfiles"] = dict(
    constants_filename="../libs/cleoconstants.hpp",  # name of file for values of physical constants
    grid_filename=str(
        path2build / "share/eurec4a1d_ddimlessGBxboundaries.dat"
    ),  # binary filename for initialisation of GBxs / GbxMaps
)

updated_configfile["initsupers"] = dict(
    type="frombinary",
    initsupers_filename=str(path2build / "share/eurec4a1d_dimlessSDsinit.dat"),
    initnsupers=0,  # Modify later!!!!
)

gridfile = updated_configfile["inputfiles"]["grid_filename"]
initSDsfile = updated_configfile["initsupers"]["initsupers_filename"]
thermofile = updated_configfile["coupled_dynamics"]["thermo"]


# %%
### ---------------------------------------------------------------- ###
### --- SETTING UP THERMODYNAMICS AND SUPERDROPLET INITIAL SETUP --- ###
### ---------------------------------------------------------------- ###


### --- settings for 1-D gridbox boundaries --- ###
cloud_altitude = cloud_observation_config["cloud"]["altitude"][0]
# only use integer precision
cloud_altitude = int(cloud_altitude)

cloud_bottom = cloud_altitude - 100
cloud_top = cloud_altitude + 100
vertical_resolution = 100

zgrid = [
    0,
    cloud_top,
    vertical_resolution,
]  # evenly spaced zhalf coords [zmin, zmax, zdelta] [m]
xgrid = np.array([0, 20])  # array of xhalf coords [m]
ygrid = np.array([0, 20])  # array of yhalf coords [m]

air_temperature_params = cloud_observation_config["thermodynamics"]["air_temperature"][
    "parameters"
]
specific_humidity_params = cloud_observation_config["thermodynamics"][
    "specific_humidity"
]["parameters"]

### --- settings for 1-D Thermodynamics --- ###
PRESS0 = 101315  # [Pa]
TEMP0 = air_temperature_params["f_0"][0]  # [K]
TEMPlapses = (
    np.array(air_temperature_params["slopes"]) * -1e3
)  # -1e3 due to conversion from dT/dz [K/m] to -dT/dz [K/km]
qvap0 = specific_humidity_params["f_0"][0]  # [Kg/Kg]
qvaplapses = (
    np.array(specific_humidity_params["slopes"]) * -1e6
)  # -1e6 due to conversion from dvap/dz [kg/kg m^-1] to -dvap/dz [g/Kg km^-1]
qcond = 0.0  # [Kg/Kg]
WVEL = 0.0  # [m/s]
Wlength = (
    0  # [m] use constant W (Wlength=0.0), or sinusoidal 1-D profile below cloud base
)

z_split_temp = air_temperature_params["x_split"]  # [m]
z_split_qvap = specific_humidity_params["x_split"]  # [m]

# create the base of the cloud as the mean of the two splits
Zbase = np.mean([z_split_temp, z_split_qvap])  # [m]

### --- settings for initial superdroplets --- ###
npergbx = 256  # number of superdroplets per gridbox

# initial superdroplet radii (and implicitly solute masses)
rspan = [1e-7, 1e-3]  # min and max range of radii to sample [m]

# initial superdroplet attributes
psd_params = cloud_observation_config["particle_size_distribution"]["parameters"]

# settings for initial superdroplet multiplicies with ATR and Aerosol from Lohmann et. al 2016 Fig. 5.5
geomeans = psd_params["geometric_means"]
geosigs = psd_params["geometric_sigmas"]
scalefacs = psd_params["scale_factors"]
numconc = np.sum(scalefacs)

### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###
# %%
### ---------------------------------------------------------------- ###
### ----- MODIFY THE CONFIG FILE PRIOR TO CREATING INPUT FILES ----- ###
### ---------------------------------------------------------------- ###

editconfigfile.edit_config_params(str(configfile), updated_configfile)

### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###

# %%


### ---------------------------------------------------------------- ###
### ------------------- BINARY FILES GENERATION--------------------- ###
### ---------------------------------------------------------------- ###
### --- ensure build, share and bin directories exist --- ###
if path2CLEO == path2build:
    raise ValueError("build directory cannot be CLEO")
else:
    Path(path2build).mkdir(exist_ok=True)
    Path(sharepath).mkdir(exist_ok=True)
    Path(binpath).mkdir(exist_ok=True)
os.system("rm " + gridfile)
os.system("rm " + initSDsfile)
os.system("rm " + thermofile[:-4] + "*")


### ----- write gridbox boundaries binary ----- ###
cgrid.write_gridboxboundaries_binary(gridfile, zgrid, xgrid, ygrid, constsfile)
with Capturing() as grid_info:
    rgrid.print_domain_info(constsfile, gridfile)

for line in grid_info:
    print(line)
    if "domain no. gridboxes:" in line:
        grid_dimensions = np.array(
            line.split(":")[-1].replace(" ", "").split("x"), dtype=int
        )
        total_number_gridboxeds = int(np.prod(grid_dimensions))
        break

# get the total number of gridboxes
updated_configfile["domain"] = dict(
    nspacedims=1,
    ngbxs=total_number_gridboxeds,
)

# update the config file
editconfigfile.edit_config_params(str(configfile), updated_configfile)

# %%

### ----- write thermodynamics binaries ----- ###
thermodyngen = thermogen.ConstHydrostaticLapseRates(
    configfile,
    constsfile,
    PRESS0,
    TEMP0,
    qvap0,
    Zbase,
    TEMPlapses,
    qvaplapses,
    qcond,
    WVEL,
    None,
    None,
    Wlength,
)
cthermo.write_thermodynamics_binary(
    thermofile, thermodyngen, configfile, constsfile, gridfile
)

### ----- write initial superdroplets binary ----- ###
nsupers = crdgens.nsupers_at_domain_top(gridfile, constsfile, npergbx, cloud_bottom)

# get total number of superdroplets
total_nsupers = int(np.sum(list(nsupers.values())))
updated_configfile["initsupers"]["initnsupers"] = total_nsupers
# update the config file
editconfigfile.edit_config_params(str(configfile), updated_configfile)


# create initial superdroplets coordinates
coord3gen = crdgens.SampleCoordGen(True)  # sample coord3 randomly
coord1gen = None  # do not generate superdroplet coord2s
coord2gen = None  # do not generate superdroplet coord2s

# create initial superdroplets attributes
xiprobdist = probdists.LnNormal(geomeans, geosigs, scalefacs)
radiigen = rgens.SampleLog10RadiiGen(rspan)

# create uniform dry radii
monodryr = 1e-9  # all SDs have this same dryradius [m]
dryradiigen = rgens.MonoAttrGen(monodryr)

# write initial superdroplets binary
initattrsgen = attrsgen.AttrsGenerator(
    radiigen, dryradiigen, xiprobdist, coord3gen, coord1gen, coord2gen
)
csupers.write_initsuperdrops_binary(
    initSDsfile, initattrsgen, configfile, constsfile, gridfile, nsupers, numconc
)
### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###

### ---------------------------------------------------------------- ###
### --------- MODIFY THE CONFIG FILE PRIOR TO RUNNING CLEO --------- ###
### ---------------------------------------------------------------- ###

editconfigfile.edit_config_params(str(configfile), updated_configfile)

### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###
### ---------------------- PLOT INIT FIGURES ----------------------- ###
### ---------------------------------------------------------------- ###

isfigures = [True, True]  # booleans for [making, saving] initialisation figures
savefigpath = (
    str(rawdirectory_individual / "figures") + "/"
)  # directory for saving figures

SDgbxs2plt = list(range(total_number_gridboxeds - 2, total_number_gridboxeds - 1))
SDgbxs2plt = [random.choice(SDgbxs2plt)]  # choose random gbx from list to plot

### ----- show (and save) plots of binary file data ----- ###
if isfigures[0]:
    if isfigures[1]:
        Path(savefigpath).mkdir(exist_ok=True)
    rgrid.plot_gridboxboundaries(constsfile, gridfile, savefigpath, isfigures[1])
    rthermo.plot_thermodynamics(
        constsfile, configfile, gridfile, thermofile, savefigpath, isfigures[1]
    )
    rsupers.plot_initGBxs_distribs(
        configfile,
        constsfile,
        initSDsfile,
        gridfile,
        savefigpath,
        isfigures[1],
        SDgbxs2plt,
    )
### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###


### ---------------------------------------------------------------- ###
### ---------------------- RUN CLEO EXECUTABLE --------------------- ###
### ---------------------------------------------------------------- ###
os.chdir(path2build)
os.system("pwd")
os.system("rm -rf " + dataset)  # delete any existing dataset
executable = str(path2build) + "/examples/eurec4a1d/src/eurec4a1D"
print("Executable: " + executable)
print("Config file: " + str(configfile))
os.system(executable + " " + str(configfile))
### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###


# %%
### ---------------------------------------------------------------- ###
### ---------------- CONVERT OUTPUT TO REGULAR ARRAY --------------- ###
### ---------------------------------------------------------------- ###
OUTPUT_FILEPATH = rawdirectory_individual / "full_dataset.nc"
# OUTPUT_FILEPATH = '/home/m/m301096/CLEO/data/output/raw/no_aerosols/cluster_18/clusters_18/full_dataset.nc'
### ----------------------- INPUT PARAMETERS ----------------------- ###
### --- essential paths and filenames --- ###

# read in constants and intial setup from setup .txt file
config = pysetuptxt.get_config(setupfile, nattrs=3, isprint=False)
consts = pysetuptxt.get_consts(setupfile, isprint=False)
# Create a first simple dataset to have the coordinates for later netcdf creation
sddata = pyzarr.get_supers(str(dataset), consts)
lagrange = sddata.to_Dataset(check_indices_uniqueness=True)
lagrange.to_netcdf(OUTPUT_FILEPATH)

### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###

# %%
