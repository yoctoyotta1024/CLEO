"""
----- CLEO -----
File: create_gbxboundariesbinary_script.py
Project: scripts
Created Date: Tuesday 24th October 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
Last Modified: Monday 17th June 2024
Modified By: CB
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
from pySD.gbxboundariesbinary_src.create_gbxboundaries import (
    write_gridboxboundaries_binary,
)
from pySD.gbxboundariesbinary_src.read_gbxboundaries import (
    print_domain_info,
    plot_gridboxboundaries,
)

### ----------------------- INPUT PARAMETERS ----------------------- ###
### absolute or relative paths for build and CLEO directories
path2CLEO = sys.argv[1]
path2build = sys.argv[2]
configfile = sys.argv[3]

# booleans for [making, saving] initialisation figures
isfigures = [True, True]

### essential paths and filenames
constsfile = path2CLEO + "/libs/cleoconstants.hpp"
binariespath = path2build + "/share/"
savefigpath = path2build + "/bin/"

gridfile = (
    binariespath + "/dimlessGBxboundaries.dat"
)  # note this should match config.yaml

### input parameters for zcoords of gridbox boundaries
zmax = 20  # maximum z coord [m]
zmin = 0  # minimum z coord [m]
zdelta = 20  # even spacing
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
    Path(path2build).mkdir(exist_ok=True)
    Path(binariespath).mkdir(exist_ok=True)

### write gridbox boundaries binary
write_gridboxboundaries_binary(gridfile, zgrid, xgrid, ygrid, constsfile)
print_domain_info(constsfile, gridfile)

### plot gridbox boundaries binary
if isfigures[0]:
    if isfigures[1]:
        Path(savefigpath).mkdir(exist_ok=True)
    plot_gridboxboundaries(constsfile, gridfile, savefigpath, isfigures[1])
### ---------------------------------------------------------------- ###
