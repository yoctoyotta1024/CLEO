"""
----- CLEO -----
File: create_initsuperdropsbinary_script.py
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
example of various ways to use cleopy module to create binary file for
initial superdroplet conditions to read into CLEO SDM
"""

import argparse
from pathlib import Path
# import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument(
    "path2CLEO", type=Path, help="Absolute path to CLEO directory (for cleopy)"
)
parser.add_argument("path2build", type=Path, help="Absolute path to build directory")
parser.add_argument(
    "config_filename", type=Path, help="Absolute path to configuration YAML file"
)
args = parser.parse_args()

from cleopy import geninitconds
from cleopy.initsuperdropsbinary_src import rgens, probdists, attrsgen, crdgens

### ----------------------- INPUT PARAMETERS ----------------------- ###
### --- absolute or relative paths for --- ###
### ---   build and CLEO directories --- ###
path2CLEO = args.path2CLEO
path2build = args.path2build
config_filename = args.config_filename

# booleans for [showing, saving] initialisation figures
isfigures = [True, True]
gbxs2plt = "all"  # indexes of GBx index of SDs to plot (nb. "all" can be very slow)

### essential paths and filenames
constants_filename = path2CLEO / "libs" / "cleoconstants.hpp"
binariespath = path2build / "share"
savefigpath = path2build / "bin"

grid_filename = (
    binariespath / "dimlessGBxboundaries.dat"
)  # note this should match config.yaml
initsupers_filename = (
    binariespath / "dimlessSDsinit.dat"
)  # note this should match config.yaml

### --- Number of Superdroplets per Gridbox --- ###
### ---        (an int or dict of ints)     --- ###
# zlim = 800
# npergbx = 8192
# nsupers =  crdgens.nsupers_at_domain_base(grid_filename, constants_filename, npergbx, zlim) # supers where z <= zlim
# nsupers = crdgens.nsupers_at_domain_top(
#     grid_filename, constants_filename, npergbx, zlim
# )  # supers where z >= zlim
nsupers = 256
### ------------------------------------------- ###

### --- Choice of Superdroplet Radii Generator --- ###
# monor                = 0.05e-6                        # all SDs have this same radius [m]
# radiigen  =  rgens.MonoAttrGen(monor)                 # all SDs have the same radius [m]

rspan = [1e-6, 1e-3]  # min and max range of radii to sample [m]
radiigen = rgens.SampleLog10RadiiGen(rspan)  # radii are sampled from rspan [m]
### ---------------------------------------------- ###

### --- Choice of Superdroplet Dry Radii Generator --- ###
monodryr = 1e-9  # all SDs have this same dryradius [m]
dryradiigen = rgens.MonoAttrGen(monodryr)  # all SDs have the same dryradius [m]

# dryr_sf = 1.0  # scale factor for dry radii [m]
# dryradiigen = dryrgens.ScaledRadiiGen(dryr_sf)  # dryradii are 1/sf of radii [m]

### ---------------------------------------------- ###

### --- Choice of Droplet Radius Probability Distribution --- ###
# dirac0               = monor                         # radius in sample closest to this value is dirac delta peak
# numconc              = 1e6                         # total no. conc of real droplets [m^-3]
# numconc              = 512e6                         # total no. conc of real droplets [m^-3]
# xiprobdist = probdists.DiracDelta(dirac0)

# geomeans           = [0.075e-6]                  # lnnormal modes' geometric mean droplet radius [m]
# geosigs            = [1.5]                       # lnnormal modes' geometric standard deviation
# scalefacs          = [1]                         # relative heights of modes
# geomeans             = [0.02e-6, 0.2e-6, 3.5e-6]
# geosigs              = [1.55, 2.3, 2]
# scalefacs            = [1, 0.3, 0.025]
# geomeans = [0.02e-6, 0.15e-6]
# geosigs = [1.4, 1.6]
# scalefacs = [0.6, 0.4]
# numconc = np.sum(scalefacs) * 1e9
# xiprobdist = probdists.LnNormal(geomeans, geosigs, scalefacs)

# volexpr0             = 30.531e-6                   # peak of volume exponential distribution [m]
# numconc              = 2**(23)                     # total no. conc of real droplets [m^-3]
# xiprobdist = probdists.VolExponential(volexpr0, rspan)

reff = 7e-6  # effective radius [m]
nueff = 0.08  # effective variance
# xiprobdist = probdists.ClouddropsHansenGamma(reff, nueff)
rdist1 = probdists.ClouddropsHansenGamma(reff, nueff)
nrain = 3000  # raindrop concentration [m^-3]
qrain = 0.9  # rainwater content [g/m^3]
dvol = 8e-4  # mean volume diameter [m]
# xiprobdist = probdists.RaindropsGeoffroyGamma(nrain, qrain, dvol)
rdist2 = probdists.RaindropsGeoffroyGamma(nrain, qrain, dvol)
numconc = 1e9  # [m^3]
distribs = [rdist1, rdist2]
scalefacs = [1000, 1]
xiprobdist = probdists.CombinedRadiiProbDistribs(distribs, scalefacs)

### --------------------------------------------------------- ###

### --- Choice of Superdroplet Coord3 Generator --- ###
# monocoord3           = 1000                        # all SDs have this same coord3 [m]
# coord3gen            =  crdgens.MonoCoordGen(monocoord3)
coord3gen = crdgens.SampleCoordGen(True)  # sample coord3 range randomly or not
# coord3gen = None  # do not generate superdroplet coord3s
### ----------------------------------------------- ###

### --- Choice of Superdroplet Coord1 Generator --- ###
# monocoord1           = 200                        # all SDs have this same coord1 [m]
# coord1gen            =  crdgens.MonoCoordGen(monocoord1)
# coord1gen = crdgens.SampleCoordGen(True)  # sample coord1 range randomly or not
coord1gen = None  # do not generate superdroplet coord1s
### ----------------------------------------------- ###

### --- Choice of Superdroplet Coord2 Generator --- ###
# monocoord2           = 1000                        # all SDs have this same coord2 [m]
# coord2gen            =  crdgens.MonoCoordGen(monocoord2)
# coord2gen = crdgens.SampleCoordGen(True)  # sample coord1 range randomly or not
coord2gen = None  # do not generate superdroplet coord2s
### ----------------------------------------------- ###

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


### write initial superdrops binary
initattrsgen = attrsgen.AttrsGenerator(
    radiigen, dryradiigen, xiprobdist, coord3gen, coord1gen, coord2gen
)
geninitconds.generate_initial_superdroplet_conditions(
    initattrsgen,
    initsupers_filename,
    config_filename,
    constants_filename,
    grid_filename,
    nsupers,
    numconc,
    isprintinfo=True,
    isfigures=isfigures,
    savefigpath=savefigpath,
    gbxs2plt=gbxs2plt,
)
### ---------------------------------------------------------------- ###
