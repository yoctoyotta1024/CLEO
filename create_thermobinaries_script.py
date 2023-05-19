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

### booleans for [making+showing, saving] figures
isfigures = [True, True]

### essential paths and filenames
constsfile = path2CLEO+"libs/claras_SDconstants.hpp"
configfile = path2CLEO+"src/config/config.txt"
binariespath = path2build+"/share/"
savefigpath = path2build+"/bin/"

gridfile =  binariespath+"/dimlessGBxboundaries.dat" # note this should match config.txt
thermofile =  binariespath+"dimlessthermodynamics.dat"


### --- Choose Initial Thermodynamic Conditions for Gridboxes  --- ###

# ### --- Constant and Uniform --- ###
# P_INIT = 100000.0                       # initial pressure [Pa]
# TEMP_INIT = 273.15                      # initial parcel temperature [T]
# relh_init = 95.0                        # initial relative humidity (%)
# qc_init = 0.0                           # initial liquid water content []
# W_INIT = 0.0                            # initial vertical (z) velocity [m/s]
# U_INIT = 0.0                            # initial horizontal x velocity [m/s]
# V_INIT = 0.0                            # initial horizontal y velocity [m/s]
# gen = thermogen.ConstUniformThermo(P_INIT, TEMP_INIT, None,
#                                     qc_init, W_INIT, U_INIT, V_INIT,
#                                     relh=relh_init, constsfile=constsfile)

### --- 2D Flow Field with Hydrostatic --- ###
### ---       or Simple z Profile      --- ###
PRESS0 = 101500 # [Pa]
THETA = 289 # [K]
qcond = 0.0 # [Kg/Kg]
WMAX = 0.6 # [m/s]
VVEL = None # [m/s]
Zlength = 1500 # [m]
Xlength = 1500 # [m]

# qvapmethod = "sratio"
qvapmethod = "qvap"
sratios = [0.85, 1.05]
sratios = [0.0075, 0.0075]
zbase = 750

# qvap = 0.0075 # [Kg/Kg]
# qvapmethod = "sratio"
gen = thermogen.ConstHydrostaticAdiabat(configfile, constsfile, PRESS0, 
                                        THETA, qvapmethod, sratios, zbase,
                                        qcond, WMAX, Zlength, Xlength,
                                        VVEL)
# gen = thermogen.SimpleThermo2Dflowfield(configfile, constsfile, PRESS0,
#                                         THETA, qvapmethod, sratios, zbase,
#                                         qcond, WMAX, Zlength, Xlength,
#                                         VVEL)
### ---------------------------------------------------------------- ###

### -------------------- BINARY FILE GENERATION--------------------- ###
cthermo.write_thermodynamics_binary(thermofile, gen, configfile,
                                    constsfile, gridfile)

if isfigures[0]:
    if isfigures[1]:
        Path(savefigpath).mkdir(exist_ok=True) 
    rthermo.plot_thermodynamics(constsfile, configfile, gridfile,
                                          thermofile, savefigpath,
                                          isfigures[1])
### ---------------------------------------------------------------- ###