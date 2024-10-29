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

import os
import shutil
import subprocess
import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

path2CLEO = Path(sys.argv[1])
path2build = Path(sys.argv[2])
configfile = Path(sys.argv[3])
outputdir = Path(sys.argv[4])
buildtype = sys.argv[5]
nruns = 2

sys.path.append(str(path2CLEO))  # imports from pySD
sys.path.append(
    str(path2CLEO / "examples" / "exampleplotting")
)  # imports from example plots package


from pySD.gbxboundariesbinary_src import read_gbxboundaries as rgrid
from pySD.gbxboundariesbinary_src import create_gbxboundaries as cgrid
from pySD.initsuperdropsbinary_src import crdgens, rgens, dryrgens, probdists, attrsgen
from pySD.initsuperdropsbinary_src import create_initsuperdrops as csupers
from pySD.initsuperdropsbinary_src import read_initsuperdrops as rsupers
from pySD.thermobinary_src import thermogen
from pySD.thermobinary_src import create_thermodynamics as cthermo
from pySD.thermobinary_src import read_thermodynamics as rthermo

### ---------------------------------------------------------------- ###
### ----------------------- INPUT PARAMETERS ----------------------- ###
### ---------------------------------------------------------------- ###
### --- essential paths and filenames --- ###
# path and filenames for creating initial SD conditions
constsfile = path2CLEO / "libs" / "cleoconstants.hpp"
binpath = path2build / "bin"
sharepath = path2build / "share"
gridfile = sharepath / "spd_dimlessGBxboundaries.dat"
initSDsfile = sharepath / "spd_dimlessSDsinit.dat"
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
def read_statsfile(statsfile):
    stats = {}
    with open(statsfile, "r") as file:
        for line in file:
            # Check if the line starts with '###'
            if not line.startswith("###"):
                # Process the line
                line = line.strip().split()
                stats[line[0]] = float(line[1])

    return stats


def write_outstats(nruns, n, outdatafile, buildtype, stats):
    """if outdatafile doesn't already exist, creates new file with
    a header. else appends to end of file"""

    try:
        # Try to open the file for exclusive creation
        with open(outdatafile, "x") as file:
            # Perform operations on the new file if needed
            header = "### Wall Clock time For Timestepping\n"
            header += "### columns are: "
            header += "test_run gpus_cpus/s cpus/s serial/s"
            file.write(header)
        print(f"--- new stats output: '{outdatafile}' created ---")
    except FileExistsError:
        print(f"stats output file '{outdatafile}' already exists")

    # write new line at number = existing number of (non-header) lines + 1
    with open(outdatafile, "r") as file:
        lines = file.readlines()

    if buildtype == "cuda":
        line = "\n" + str(n) + " " + str(stats["tstep"])
        with open(outdatafile, "a") as file:
            file.write(line)

    else:
        nline = len(lines) - nruns + n
        lines[nline] = lines[nline].rstrip() + " " + str(stats["tstep"]) + "\n"
        with open(outdatafile, "w") as file:
            file.writelines(lines)


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
shutil.rmtree(gridfile, ignore_errors=True)
shutil.rmtree(initSDsfile, ignore_errors=True)
all_thermofiles = thermofiles.parent / Path(f"{thermofiles.stem}*{thermofiles.suffix}")
shutil.rmtree(all_thermofiles, ignore_errors=True)

### ----- write gridbox boundaries binary ----- ###
cgrid.write_gridboxboundaries_binary(gridfile, zgrid, xgrid, ygrid, constsfile)
rgrid.print_domain_info(constsfile, gridfile)

### ----- write thermodynamics binaries ----- ###
thermodyngen = thermogen.SimpleThermo2DFlowField(
    configfile,
    constsfile,
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
cthermo.write_thermodynamics_binary(
    thermofiles, thermodyngen, configfile, constsfile, gridfile
)


### ----- write initial superdroplets binary ----- ###
nsupers = crdgens.nsupers_at_domain_base(gridfile, constsfile, npergbx, zlim)
coord3gen = crdgens.SampleCoordGen(True)  # sample coord3 randomly
coord1gen = crdgens.SampleCoordGen(True)  # sample coord1 randomly
coord2gen = crdgens.SampleCoordGen(True)  # sample coord2 randomly
xiprobdist = probdists.LnNormal(geomeans, geosigs, scalefacs)
radiigen = rgens.SampleLog10RadiiGen(rspan)  # randomly sample radii from rspan [m]
dryradiigen = dryrgens.ScaledRadiiGen(1.0)

initattrsgen = attrsgen.AttrsGenerator(
    radiigen, dryradiigen, xiprobdist, coord3gen, coord1gen, coord2gen
)
csupers.write_initsuperdrops_binary(
    initSDsfile, initattrsgen, configfile, constsfile, gridfile, nsupers, numconc
)

### ----- show (and save) plots of binary file data ----- ###
if isfigures[0]:
    rgrid.plot_gridboxboundaries(constsfile, gridfile, savefigpath, isfigures[1])
    rthermo.plot_thermodynamics(
        constsfile, configfile, gridfile, thermofiles, savefigpath, isfigures[1]
    )
    rsupers.plot_initGBxs_distribs(
        configfile,
        constsfile,
        initSDsfile,
        gridfile,
        savefigpath,
        isfigures[1],
        SDgbxs2plt,
    )
    plt.close()
### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###

### ---------------------------------------------------------------- ###
### ---------------------- RUN CLEO EXECUTABLE --------------------- ###
### ---------------------------------------------------------------- ###
executable = path2build / "examples" / "speedtest" / "src" / "spdtest"
for n in range(nruns):
    os.chdir(path2build)
    shutil.rmtree(dataset, ignore_errors=True)  # delete any existing dataset
    print("Executable: " + str(executable))
    print("Config file: " + str(configfile))
    subprocess.run([executable, configfile])

    # copy speed results to new file
    print("--- reading runtime statistics ---")
    stats = read_statsfile(statsfile)
    for key, value in stats.items():
        print(key + ": {:.3f}s".format(value))
    write_outstats(nruns, n, outdatafile, buildtype, stats)
    print("--- runtime stats written to file ---")
### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###
