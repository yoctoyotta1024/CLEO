"""
----- CLEO -----
File: speedtest.py
Project: speedtest
Created Date: Friday 17th November 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
Last Modified: Sunday 14th April 2024
Modified By: CB
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
Copyright (c) 2023 MPI-M, Clara Bayley
-----
File Description:
Script generates input files, then runs CLEO executable "spdtest" to check performance
of CLEO usign different build configurations (e.g. serial, OpenmP and CUDA parallelism).
"""

import glob
import os
import shutil
import subprocess
import sys
import numpy as np
from pathlib import Path

path2CLEO = Path(sys.argv[1])
path2build = Path(sys.argv[2])
config_filename = Path(sys.argv[3])
outputdir = Path(sys.argv[4])
path2kokkostools = Path(sys.argv[5])
buildtype = sys.argv[6]
nruns = 2

sys.path.append(str(path2CLEO))  # imports from pySD
sys.path.append(
    str(path2CLEO / "examples" / "exampleplotting")
)  # imports from example plots package

from pySD import geninitconds
from pySD.initsuperdropsbinary_src import crdgens, rgens, dryrgens, probdists, attrsgen
from pySD.thermobinary_src import thermogen

### ---------------------------------------------------------------- ###
### ----------------------- INPUT PARAMETERS ----------------------- ###
### ---------------------------------------------------------------- ###
### --- essential paths and filenames --- ###
# path and filenames for creating initial SD conditions
constants_filename = path2CLEO / "libs" / "cleoconstants.hpp"
binpath = path2build / "bin"
sharepath = path2build / "share"
grid_filename = sharepath / "spd_dimlessGBxboundaries.dat"
initsupers_filename = sharepath / "spd_dimlessSDsinit.dat"
thermofiles = sharepath / "spd_dimlessthermo.dat"

# path and file names for plotting results
setupfile = binpath / "spd_setup.txt"
statsfile = binpath / "spd_stats.txt"
dataset = binpath / "spd_sol.zarr"

### --- plotting initialisation figures --- ###
isfigures = [False, False]  # booleans for [making, saving] initialisation figures
savefigpath = outputdir  # directory for saving figures
SDgbxs2plt = [0]  # gbxindex of SDs to plot (nb. "all" can be very slow)
outdatafile = outputdir / "spd_allstats.txt"  # file to write out stats to

### --- settings for 3-D gridbox boundaries --- ###
zgrid = [0, 1500, 50]  # evenly spaced zhalf coords [zmin, zmax, zdelta] [m]
xgrid = [0, 1500, 50]  # evenly spaced xhalf coords [m]
ygrid = np.array([0, 25, 50])  # array of yhalf coords [m]

### --- settings for initial superdroplets --- ###
# settings for initial superdroplet coordinates
zlim = 1500  # max z coord of superdroplets
npergbx = 4  # number of superdroplets per gridbox

# [min, max] range of initial superdroplet radii (and implicitly solute masses)
rspan = [3e-9, 3e-6]  # [m]

# settings for initial superdroplet multiplicies
# (from bimodal Lognormal distribution)
geomeans = [0.02e-6, 0.15e-6]
geosigs = [1.4, 1.6]
scalefacs = [6e6, 4e6]
numconc = np.sum(scalefacs)

### --- settings for 3D Thermodynamics --- ###
PRESS0 = 100000  # [Pa]
THETA = 298.15  # [K]
qcond = 0.0  # [Kg/Kg]
WMAX = 3.0  # [m/s]
VVEL = 1.0  # [m/s]
Zlength = 1500  # [m]
Xlength = 1500  # [m]
qvapmethod = "sratio"
Zbase = 750  # [m]
sratios = [0.85, 1.1]  # s_ratio [below, above] Zbase
### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###


### ---------------------------------------------------------------- ###
### --------------------- FUNCTION DEFINITIONS --------------------- ###
### ---------------------------------------------------------------- ###
class KpKernelTimer:
    def __init__(self, kokkos_tools_lib: Path):
        self.kokkos_tools_lib = kokkos_tools_lib
        self.kp_reader = self.kokkos_tools_lib / ".." / "bin" / "kp_reader"

        os.environ["KOKKOS_TOOLS_LIBS"] = str(
            self.kokkos_tools_lib / "libkp_kernel_timer.so"
        )
        print("Using Kokkos Profiling Tool", os.environ["KOKKOS_TOOLS_LIBS"])
        print("Using Kokkos Tool Reader", self.kp_reader)

    def postprocess(self, filespath: Path, txt_filepath: Path, txt_filelabel: str):
        # Add kokkos_tools_lib to LD_LIBRARY_PATH
        ld_lib_path = os.environ.get("LD_LIBRARY_PATH", "")
        os.environ["LD_LIBRARY_PATH"] = f"{self.kokkos_tools_lib}:{ld_lib_path}"

        # Use glob to find all .dat files in the specified directory
        datfiles = glob.glob(os.path.join(filespath, "*.dat"))

        # Print the list of .dat files
        for f, filename in enumerate(datfiles):
            txt_filename = txt_filepath / Path(txt_filelabel + f"_{f}.txt")
            cmd = [str(self.kp_reader), filename]
            with open(txt_filename, "w") as wfile:
                subprocess.run(cmd, stdout=wfile, stderr=subprocess.STDOUT)


### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###

### ---------------------------------------------------------------- ###
### ------------------- BINARY FILES GENERATION--------------------- ###
### ---------------------------------------------------------------- ###
### --- ensure build, share and bin directories exist --- ###
if path2CLEO == path2build:
    raise ValueError("build directory cannot be CLEO")
else:
    path2build.mkdir(exist_ok=True)
    sharepath.mkdir(exist_ok=True)
    binpath.mkdir(exist_ok=True)
    if isfigures[1]:
        savefigpath.mkdir(exist_ok=True)

### --- delete any existing initial conditions --- ###
shutil.rmtree(grid_filename, ignore_errors=True)
shutil.rmtree(initsupers_filename, ignore_errors=True)
all_thermofiles = thermofiles.parent / Path(f"{thermofiles.stem}*{thermofiles.suffix}")
shutil.rmtree(all_thermofiles, ignore_errors=True)

### ----- write gridbox boundaries binary ----- ###
geninitconds.generate_gridbox_boundaries(
    grid_filename,
    zgrid,
    xgrid,
    ygrid,
    constants_filename,
    isfigures=isfigures,
    savefigpath=savefigpath,
)

### ----- write thermodynamics binaries ----- ###
thermodyngen = thermogen.SimpleThermo2DFlowField(
    config_filename,
    constants_filename,
    PRESS0,
    THETA,
    qvapmethod,
    sratios,
    Zbase,
    qcond,
    WMAX,
    Zlength,
    Xlength,
    VVEL,
)
geninitconds.generate_thermodynamics_conditions_fromfile(
    thermofiles,
    thermodyngen,
    config_filename,
    constants_filename,
    grid_filename,
    isfigures=isfigures,
    savefigpath=savefigpath,
)

### ----- write initial superdroplets binary ----- ###
nsupers = crdgens.nsupers_at_domain_base(
    grid_filename, constants_filename, npergbx, zlim
)
coord3gen = crdgens.SampleCoordGen(True)  # sample coord3 randomly
coord1gen = crdgens.SampleCoordGen(True)  # sample coord1 randomly
coord2gen = crdgens.SampleCoordGen(True)  # sample coord2 randomly
xiprobdist = probdists.LnNormal(geomeans, geosigs, scalefacs)
radiigen = rgens.SampleLog10RadiiGen(rspan)  # randomly sample radii from rspan [m]
dryradiigen = dryrgens.ScaledRadiiGen(1.0)

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
    isfigures=isfigures,
    savefigpath=savefigpath,
    gbxs2plt=SDgbxs2plt,
)
### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###

### ---------------------------------------------------------------- ###
### ---------------------- RUN CLEO EXECUTABLE --------------------- ###
### ---------------------------------------------------------------- ###
profiler = KpKernelTimer(path2kokkostools)
executable = path2build / "examples" / "speedtest" / "src" / "spdtest"
for n in range(nruns):
    os.chdir(path2build)
    shutil.rmtree(dataset, ignore_errors=True)  # delete any existing dataset
    print("Executable: " + str(executable))
    print("Config file: " + str(config_filename))
    subprocess.run([executable, config_filename])
profiler.postprocess(path2build, outputdir, buildtype)
### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###
