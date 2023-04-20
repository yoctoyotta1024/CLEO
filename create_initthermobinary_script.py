import numpy as np
from pySD.thermobinary_src import thermogen
from pySD.thermobinary_src import create_thermodynamics as cthermo
from pySD.thermobinary_src import read_thermodynamics as rthermo

abspath = "/Users/yoctoyotta1024/Documents/b1_springsummer2023/CLEO/"
constsfile = abspath+"libs/claras_SDconstants.hpp"
configfile = abspath+"src/config/config.txt"

spath = abspath+"build/share/"
gridfile = spath+"dimlessGBxboundaries.dat"
thermofile =  spath+"dimlessthermodynamics.dat"
binpath = abspath+"build/bin/"

### booleans for [making+showing, saving] figures
isfigures = [True, True]

### initial thermodynamic conditions for all gridboxes ###
P_INIT = 100000.0                       # initial pressure [Pa]
TEMP_INIT = 273.15                      # initial parcel temperature [T]
relh_init = 95.0                        # initial relative humidity (%)
qc_init = 0.0                           # initial liquid water content []
W_INIT = 0.0                            # initial vertical (z) velocity [m/s]
U_INIT = None                           # initial horizontal x velocity [m/s]
V_INIT = None                           # initial horizontal y velocity [m/s]

thermogen = thermogen.ConstUniformThermo(P_INIT, TEMP_INIT, relh_init,
                                       qc_init, W_INIT, U_INIT, V_INIT,
                                       constsfile)
cthermo.write_thermodynamics_binary(thermofile, thermogen, configfile,
                                    constsfile, gridfile)

if isfigures[0]:
    #t2plt = 0.0
    t2plt =  "all"
    rthermo.plot_thermodynamics_timeslice(constsfile, configfile, gridfile,
                                          thermofile, t2plt, binpath,
                                          isfigures[1])
