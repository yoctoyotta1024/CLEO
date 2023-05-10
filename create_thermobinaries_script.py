import numpy as np
from pySD.thermobinary_src import thermogen
from pySD.thermobinary_src import create_thermodynamics as cthermo
from pySD.thermobinary_src import read_thermodynamics as rthermo

abspath = "/Users/yoctoyotta1024/Documents/b1_springsummer2023/CLEO/"
# abspath = "/home/m/m300950/CLEO/"
#abspath = sys.argv[1]
constsfile = abspath+"libs/claras_SDconstants.hpp"
configfile = abspath+"src/config/config.txt"

spath = abspath+"build/share/"
gridfile = spath+"dimlessGBxboundaries.dat"
thermofile =  spath+"dimlessthermodynamics.dat"
binpath = abspath+"build/bin/"

### booleans for [making+showing, saving] figures
isfigures = [True, True]

### initial thermodynamic conditions for all gridboxes ###

# # ----- Constant and Uniform ----- #
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

# ----- 2D Flow Field with Hydrostatic or Simple z Profile ----- #
PRESS0 = 101500 # [Pa]
THETA = 289 # [K]
qcond = 0.0 # [Kg/Kg]
WMAX = 0.6 # [m/s]
VVEL = None # [m/s]
Zlength = 1500 # [m]
Xlength = 1500 # [m]

qvapmethod, sratio = "sratio", 0.85
zbase = 750
# gen = thermogen.SimpleThermo2Dflowfield(configfile, constsfile, PRESS0,
#                                         THETA, "sratio", zbase, sratio,
#                                         qcond, WMAX, Zlength, Xlength,
#                                         VVEL)

qvap = 0.0075 # [Kg/Kg]
gen = thermogen.ConstHydrostaticAdiabat(configfile, constsfile, PRESS0, 
                                        THETA, qvap, qcond, WMAX, 
                                        Zlength, Xlength, VVEL)
# -------------------------------------------------------------- #

cthermo.write_thermodynamics_binary(thermofile, gen, configfile,
                                    constsfile, gridfile)

if isfigures[0]:
    rthermo.plot_thermodynamics(constsfile, configfile, gridfile,
                                          thermofile, binpath,
                                          isfigures[1])