import sys
import numpy as np
from pathlib import Path

from pySD.gbxboundariesbinary_src.create_gbxboundaries import *
from pySD.gbxboundariesbinary_src.read_gbxboundaries import *

### path and filenames
# path2CLEO = "/Users/yoctoyotta1024/Documents/b1_springsummer2023/CLEO/"
# path2CLEO = "/home/m/m300950/CLEO/"
path2CLEO = sys.argv[1]
constsfile = path2CLEO+"libs/claras_SDconstants.hpp"
configfile = path2CLEO+"src/config/config.txt"

gridfile = path2CLEO+"/build/share/dimlessGBxboundaries.dat"
binpath = path2CLEO+"/build/bin/"

### booleans for [making+showing, saving] figures
isfigures = [True, True]

### input parameters for zcoords of gridbox boundaries
zmax = 1500 # maximum z coord [m]
zmin = 0 # minimum z coord [m]
zdelta = 50 # even spacing
zgrid = [zmin, zmax, zdelta] 

### input parameters for x coords of gridbox boundaries
xgrid = [0, 1500, 50]

### input parameters for y coords of gridbox boundaries
ygrid = np.asarray([0, 50])
# ygrid = [0, 1500, 200]

Path(gridfilepath).mkdir(exist_ok=True) 
Path(binpath).mkdir(exist_ok=True) 
write_gridboxboundaries_binary(gridfile, zgrid, xgrid, ygrid, constsfile)
print_domain_info(constsfile, gridfile)

if isfigures[0]:
  plot_gridboxboundaries(constsfile, gridfile, binpath, isfigures[1])


