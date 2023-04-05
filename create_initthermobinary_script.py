import numpy as np
from pySD.initthermobinary_src.create_initthermo import *

abspath = "/Users/yoctoyotta1024/Documents/b1_springsummer2023/CLEO/"
constsfile = abspath+"libs/claras_SDconstants.hpp"
gridfile = abspath+"build/share/dimlessGBxboundaries.dat"

### initial thermodynamic conditions for all gridboxes ###
P_INIT = 100000.0                       # initial pressure [Pa]
TEMP_INIT = 273.15                      # initial parcel temperature [T]
relh_init = 95.0                        # initial relative humidity (%)
qc_init = 0.0                           # initial liquid water content []

from pySD.gbxboundariesbinary_src import read_gbxboundaries
zhalf, xhalf, yhalf = read_gbxboundaries.get_gridboxboundaries(constsfile,
                                                               gridfile)

print(zhalf)
print(xhalf)
print(yhalf)
