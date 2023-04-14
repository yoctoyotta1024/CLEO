import sys
import numpy as np
from pathlib import Path

from pySD.gbxboundariesbinary_src.create_gbxboundaries import *
from pySD.gbxboundariesbinary_src.read_gbxboundaries import *

### path and filenames
abspath = "/Users/yoctoyotta1024/Documents/b1_springsummer2023/CLEO/"
#abspath = sys.argv[1]
constsfile = abspath+"libs/claras_SDconstants.hpp"
configfile = abspath+"src/config/config.txt"

gridfilepath = abspath+"build/share/"
gridfile = gridfilepath+"dimlessGBxboundaries.dat"
binpath = abspath+"build/bin/"

### booleans for [making+showing, saving] figures
isfigures = [True, True]

### input parameters for zcoords of gridbox boundaries
zmax = 4000 # maximum z coord [m]
zmin = 0 # minimum z coord [m]
zdelta = 1000 # even spacing
zgrid = [zmin, zmax, zdelta] 
#zgrid = np.asarray([0, 500, 1000, 2500, 5000])

### input parameters for x coords of gridbox boundaries
xgrid = np.asarray([0, 1000])

### input parameters for y coords of gridbox boundaries
ygrid = np.asarray([0, 1000])

Path(gridfilepath).mkdir(exist_ok=True) 
Path(binpath).mkdir(exist_ok=True) 
write_gridboxboundaries_binary(gridfile, zgrid, xgrid, ygrid, constsfile)
print_domain_info(constsfile, gridfile)

if isfigures[0]:
  plot_gridboxboundaries(constsfile, gridfile, binpath, isfigures[1])


