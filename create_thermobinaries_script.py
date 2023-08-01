import sys
import numpy as np
from pathlib import Path

from pySD.thermobinary_src import thermogen
from pySD.thermobinary_src import create_thermodynamics as cthermo
from pySD.thermobinary_src import read_thermodynamics as rthermo

### ----------------------- INPUT PARAMETERS ----------------------- ###
### --- absolute or relative paths for --- ###
### ---   build and CLEO directories --- ###
path2CLEO = sys.argv[1]
path2build = sys.argv[2]
configfile = sys.argv[3]

### booleans for [making+showing, saving] figures
isfigures = [True, True]

### essential paths and filenames
constsfile = path2CLEO+"libs/claras_SDconstants.hpp"
binariespath = path2build+"/share/"
savefigpath = path2build+"/bin/"

gridfile =  binariespath+"/dimlessGBxboundaries.dat" # note this should match config.txt
thermofile =  binariespath+"/dimlessthermo.dat"

### --- Choose Initial Thermodynamic Conditions for Gridboxes  --- ###

### --- Constant and Uniform --- ###
P_INIT = 100000.0                       # initial pressure [Pa]
TEMP_INIT = 273.15                      # initial parcel temperature [T]
relh_init = 95.0                        # initial relative humidity (%)
qc_init = 0.0                           # initial liquid water content []
W_INIT = 0.0                            # initial vertical (z) velocity [m/s]
U_INIT = 0.0                            # initial horizontal x velocity [m/s]
V_INIT = 0.0                            # initial horizontal y velocity [m/s]
thermodyngen = thermogen.ConstUniformThermo(P_INIT, TEMP_INIT, None,
                                    qc_init, W_INIT, U_INIT, V_INIT,
                                    relh=relh_init, constsfile=constsfile)

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
# # moistlayer = False
# moistlayer = {
#     "z1": 700,
#     "z2": 800,
#     "x1": 0,
#     "x2": 750,
#     "mlsratio": 1.005
# }
# thermodyngen = thermogen.ConstHydrostaticAdiabat(configfile, constsfile, PRESS0, 
#                                         THETA, qvapmethod, sratios, Zbase,
#                                         qcond, WMAX, Zlength, Xlength,
#                                         VVEL, moistlayer)
# thermodyngen = thermogen.SimpleThermo2Dflowfield(configfile, constsfile, PRESS0,
#                                         THETA, qvapmethod, sratios, Zbase,
#                                         qcond, WMAX, Zlength, Xlength,
#                                         VVEL)
### ---------------------------------------------------------------- ###

### -------------------- BINARY FILE GENERATION--------------------- ###
cthermo.write_thermodynamics_binary(thermofile, thermodyngen, configfile,
                                    constsfile, gridfile)

if isfigures[0]:
    if isfigures[1]:
        Path(savefigpath).mkdir(exist_ok=True) 
    rthermo.plot_thermodynamics(constsfile, configfile, gridfile,
                                          thermofile, savefigpath,
                                          isfigures[1])
### ---------------------------------------------------------------- ###