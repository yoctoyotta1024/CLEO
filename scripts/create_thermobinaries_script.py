"""
Copyright (c) 2024 MPI-M, Clara Bayley


----- CLEO -----
File: create_thermobinaries_script.py
Project: scripts
Created Date: Tuesday 24th October 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
example of various ways to use cleopy module to create binary files
for the dynamics to read into CLEO when using a from file data for coupled dynamics
"""

import argparse
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument(
    "path2CLEO", type=Path, help="Absolute path to CLEO directory (for cleopy)"
)
parser.add_argument("path2build", type=Path, help="Absolute path to build directory")
parser.add_argument(
    "config_filename", type=Path, help="Absolute path to configuration YAML file"
)
args = parser.parse_args()

from cleopy import geninitconds
from cleopy.thermobinary_src import thermogen, windsgen, thermodyngen

### ----------------------- INPUT PARAMETERS ----------------------- ###
### --- absolute or relative paths for --- ###
### ---   build and CLEO directories --- ###
path2CLEO = args.path2CLEO
path2build = args.path2build
config_filename = args.config_filename

# booleans for [showing, saving] initialisation figures
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

### --- Thermo (temp, press, qvap and cond) Conditions  --- ###

# ### --- Constant and Uniform --- ###
# P_INIT = 101500.0                       # initial pressure [Pa]
# TEMP_INIT = 288.15                      # initial parcel temperature [T]
# relh_init = 0.999                       # initial relative humidity (%)
# qvap = None                             # use relative humidity to set qvap
# qc_init = 0.0                           # initial liquid water content []
# thermog = thermogen.ConstUniformThermo(P_INIT, TEMP_INIT, None,
#                              qc_init, relh=relh_init,
#                              constants_filename=constants_filename)

### --- 1-D T and qv set by Lapse Rates --- ###
PRESS0 = 101315  # [Pa]
TEMP0 = 297.9  # [K]
qvap0 = 0.016  # [Kg/Kg]
Zbase = 800  # [m]
TEMPlapses = [9.8, 6.5]  # -dT/dz [K/km]
qvaplapses = [2.97, "saturated"]  # -dvap/dz [g/Kg km^-1]
qcond = 0.0  # [Kg/Kg]
thermog = thermogen.HydrostaticLapseRates(
    config_filename,
    constants_filename,
    PRESS0,
    TEMP0,
    qvap0,
    Zbase,
    TEMPlapses,
    qvaplapses,
    qcond,
)

# ### --- Hydrostatic Dry Adiabat --- ###
# ### ---   or Simple z Profile   --- ###
# PRESSz0 = 101500  # [Pa]
# THETA = 289  # [K]
# qcond = 0.0  # [Kg/Kg]
# qvapmethod = "sratio"
# Zbase = 750  # [m]
# sratios = [0.85, 1.0001]  # s_ratio [below, above] Zbase
# # moistlayer = False
# moistlayer = {"z1": 700, "z2": 800, "x1": 0, "x2": 750, "mlsratio": 1.005}
# thermog = thermogen.DryHydrostaticAdiabatic2TierRelH(
#     config_filename,
#     constants_filename,
#     PRESSz0,
#     THETA,
#     qvapmethod,
#     sratios,
#     Zbase,
#     qcond,
#     moistlayer,
# )
# thermog = thermogen.Simple2TierRelativeHumidity(config_filename, constants_filename, PRESSz0,
#                                         THETA, qvapmethod, sratios, Zbase,
#                                         qcond)

### --- Wind Field (wvel, vvel, uvel) Conditions  --- ###

### --- Constant and Uniform --- ###
W_INIT = 0.0  # initial vertical (coord3) velocity [m/s]
U_INIT = 0.0  # initial eastwards (coord1) velocity [m/s]
V_INIT = 0.0  # initial northwards (coord2) velocity [m/s]
windsg = windsgen.ConstUniformWinds(W_INIT, U_INIT, V_INIT)

# ### --- 1D Vertical Sinusoid --- ###
# WMAX = 0.0  # [m/s]
# UVEL = None
# VVEL = None
# Wlength = (
#     1000  # [m] use constant W (Wlength=0.0), or sinusoidal 1-D profile below cloud base
# )
# windsg = windsgen.SinusoidalUpdraught(WMAX, UVEL, VVEL, Wlength)

# ### --- 2D Flow Field --- ###
# WMAX = 0.6  # [m/s]
# VVEL = 1.0  # [m/s]
# Zlength = 1500  # [m]
# Xlength = 1500  # [m]
# # windsg = thermog.create_default_windsgen(
# #     WMAX, Zlength, Xlength, VVEL
# # )  # only for DryHydrostaticAdiabatic2TierRelH
# windsg = windsgen.Simple2DFlowField(
#     config_filename, constants_filename, WMAX, Zlength, Xlength, VVEL
# )

### --- Thermodynamic + Winds Conditions  --- ###
thermodyngen = thermodyngen.ThermodynamicsGenerator(thermog, windsg)
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
