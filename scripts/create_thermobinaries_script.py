"""
Copyright (c) 2024 MPI-M, Clara Bayley


----- CLEO -----
File: create_thermobinaries_script.py
Project: scripts
Created Date: Tuesday 24th October 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
Last Modified: Tuesday 7th May 2024
Modified By: CB
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
uses pySD module to create binary files
for the dynamics to read into CLEO when
using a from file data for coupled dynamics
"""

import sys
from pathlib import Path

sys.path.append(sys.argv[1])  # path to pySD (same as to CLEO)
from pySD import geninitconds
from pySD.thermobinary_src import thermogen

### ----------------------- INPUT PARAMETERS ----------------------- ###
### --- absolute or relative paths for --- ###
### ---   build and CLEO directories --- ###
path2CLEO = Path(sys.argv[1])
path2build = Path(sys.argv[2])
config_filename = Path(sys.argv[3])

# booleans for [making, saving] initialisation figures
isfigures = [True, True]

### essential paths and filenames
constants_filename = path2CLEO / "libs" / "cleoconstants.hpp"
binariespath = path2build / "share"
savefigpath = path2build / "bin"

grid_filename = (
    binariespath / "dimlessGBxboundaries.dat"
)  # note this should match config.yaml
thermofiles = binariespath / "dimlessthermo.dat"

### --- Choose Initial Thermodynamic Conditions for Gridboxes  --- ###

### --- Constant and Uniform --- ###
# P_INIT = 101500.0                       # initial pressure [Pa]
# TEMP_INIT = 288.15                      # initial parcel temperature [T]
# relh_init = 0.999                       # initial relative humidity (%)
# qc_init = 0.0                           # initial liquid water content []
# W_INIT = 0.0                           # initial vertical (coord3) velocity [m/s]
# U_INIT = 0.0                            # initial eastwards (coord1) velocity [m/s]
# V_INIT = 0.0                            # initial northwards (coord2) velocity [m/s]
# thermodyngen = thermogen.ConstUniformThermo(P_INIT, TEMP_INIT, None,
#                                     qc_init, W_INIT, U_INIT, V_INIT,
#                                     relh=relh_init, constants_filename=constants_filename)

### --- 1-D T and qv set by Lapse Rates --- ###
PRESS0 = 101315  # [Pa]
TEMP0 = 297.9  # [K]
qvap0 = 0.016  # [Kg/Kg]
Zbase = 800  # [m]
TEMPlapses = [9.8, 6.5]  # -dT/dz [K/km]
qvaplapses = [2.97, "saturated"]  # -dvap/dz [g/Kg km^-1]
qcond = 0.0  # [Kg/Kg]
WMAX = 0.0  # [m/s]
Wlength = (
    1000  # [m] use constant W (Wlength=0.0), or sinusoidal 1-D profile below cloud base
)

thermodyngen = thermogen.ConstHydrostaticLapseRates(
    config_filename,
    constants_filename,
    PRESS0,
    TEMP0,
    qvap0,
    Zbase,
    TEMPlapses,
    qvaplapses,
    qcond,
    WMAX,
    None,
    None,
    Wlength,
)

# ### --- 2D Flow Field with Hydrostatic --- ###
# ### ---       or Simple z Profile      --- ###
# PRESS0 = 101500 # [Pa]
# THETA = 289 # [K]
# qcond = 0.0 # [Kg/Kg]
# WMAX = 0.6 # [m/s]
# VVEL = None # [m/s]
# Zlength = 1500 # [m]
# Xlength = 1500 # [m]

# qvapmethod = "sratio"
# Zbase = 750 # [m]
# sratios = [0.85, 1.0001] # s_ratio [below, above] Zbase
# moistlayer = False
# moistlayer = {
#     "z1": 700,
#     "z2": 800,
#     "x1": 0,
#     "x2": 750,
#     "mlsratio": 1.005
# }
# thermodyngen = thermogen.ConstDryHydrostaticAdiabat(config_filename, constants_filename, PRESS0,
#                                         THETA, qvapmethod, sratios, Zbase,
#                                         qcond, WMAX, Zlength, Xlength,
#                                         VVEL, moistlayer)
# thermodyngen = thermogen.SimpleThermo2DFlowField(config_filename, constants_filename, PRESS0,
#                                         THETA, qvapmethod, sratios, Zbase,
#                                         qcond, WMAX, Zlength, Xlength,
#                                         VVEL)
### ---------------------------------------------------------------- ###

### -------------------- BINARY FILE GENERATION--------------------- ###
### ensure build, share and bin directories exist
if path2CLEO == path2build:
    raise ValueError("build directory cannot be CLEO")
else:
    path2build.mkdir(exist_ok=True)
    binariespath.mkdir(exist_ok=True)
    if isfigures[1]:
        savefigpath.mkdir(exist_ok=True)

geninitconds.generate_thermodynamics_conditions_fromfile(
    thermofiles,
    thermodyngen,
    config_filename,
    constants_filename,
    grid_filename,
    isfigures=isfigures,
    savefigpath=savefigpath,
)
### ---------------------------------------------------------------- ###
