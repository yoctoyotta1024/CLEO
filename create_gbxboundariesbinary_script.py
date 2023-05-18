import sys
import numpy as np
from pathlib import Path

from pySD.gbxboundariesbinary_src.create_gbxboundaries import *
from pySD.gbxboundariesbinary_src.read_gbxboundaries import *

### ----------------------- INPUT PARAMETERS ----------------------- ###
### absolute or relative paths for build and CLEO directories
path2CLEO = sys.argv[1]
path2build = sys.argv[2]

### booleans for [making+showing, saving] figures
isfigures = [True, True]

### essential paths and filenames
constsfile = path2CLEO+"libs/claras_SDconstants.hpp"
configfile = path2CLEO+"src/config/config.txt"
gridfilepath = path2build+"/share/"
savefigpath = path2build+"/bin/"

gridfile =  gridfilepath+"/dimlessGBxboundaries.dat" # note this should match config.txt

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
### ---------------------------------------------------------------- ###


### -------------------- BINARY FILE GENERATION--------------------- ###
### ensure build, share and bin directories exist
if path2CLEO == path2build:
  raise ValueError("build directory cannot be CLEO")
else:
  Path(path2build).mkdir(exist_ok=True) 
  Path(gridfilepath).mkdir(exist_ok=True) 
  Path(savefigpath).mkdir(exist_ok=True) 

### write gridbox boundaries binary
write_gridboxboundaries_binary(gridfile, zgrid, xgrid, ygrid, constsfile)
print_domain_info(constsfile, gridfile)

### plot gridbox boundaries binary
if isfigures[0]:
  plot_gridboxboundaries(constsfile, gridfile, savefigpath, isfigures[1])
### ---------------------------------------------------------------- ###