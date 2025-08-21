"""
----- CLEO -----
File: create_gbxboundariesbinary_script.py
Project: scripts
Created Date: Tuesday 24th October 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
Copyright (c) 2023 MPI-M, Clara Bayley
-----
File Description:
uses pySD module to create gridbox boundaries
binary file for input to CLEO SDM
"""

import sys
import numpy as np
from pathlib import Path

sys.path.append(sys.argv[1])  # path to pySD (same as to CLEO)
from pySD import geninitconds

### ----------------------- INPUT PARAMETERS ----------------------- ###
### absolute or relative paths for build and CLEO directories
path2CLEO = Path(sys.argv[1])
path2build = Path(sys.argv[2])
config_filename = Path(sys.argv[3])

# booleans for [showing, saving] initialisation figures
isfigures = [True, True]

### essential paths and filenames
constants_filename = path2CLEO / "libs" / "cleoconstants.hpp"
binariespath = path2build / "share"
savefigpath = path2build / "bin"

grid_filename = (
    binariespath / "dimlessGBxboundaries.dat"
)  # note this should match config.yaml

### input parameters for zcoords of gridbox boundaries
zmax = 100  # maximum z coord [m]
zmin = 0  # minimum z coord [m]
zdelta = 10  # even spacing
zgrid = [zmin, zmax, zdelta]
# zgrid = np.arange(zmin, zmax+zdelta, zdelta)

### input parameters for x coords of gridbox boundaries
xgrid = [0, 20, 20]

### input parameters for y coords of gridbox boundaries
ygrid = np.asarray([0, 20])
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
geninitconds.generate_gridbox_boundaries(
    grid_filename,
    zgrid,
    xgrid,
    ygrid,
    constants_filename,
    isprintinfo=True,
    isfigures=isfigures,
    savefigpath=savefigpath,
)
### ---------------------------------------------------------------- ###
