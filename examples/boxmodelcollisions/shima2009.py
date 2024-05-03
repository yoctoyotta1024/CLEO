'''
----- CLEO -----
File: shima2009.py
Project: boxmodelcollisions
Created Date: Friday 17th November 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
Last Modified: Friday 3rd May 2024
Modified By: CB
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
Copyright (c) 2023 MPI-M, Clara Bayley
-----
File Description:
Script generates input files, runs CLEO 0-D box model executables for
collisions with selected collision kernels (e.g. Golovin's or Long's) to create data.
Then plots results similar to Shima et al. 2009 Fig. 2
'''

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

path2CLEO = sys.argv[1]
path2build = sys.argv[2]
configfile = sys.argv[3]
kernels = sys.argv[4:]

sys.path.append(path2CLEO)  # for imports from pySD package
sys.path.append(path2CLEO+"/examples/exampleplotting/") # for imports from example plotting package

from plotssrc import shima2009fig
from pySD.sdmout_src import *
from pySD.initsuperdropsbinary_src import *
from pySD.initsuperdropsbinary_src import create_initsuperdrops as csupers
from pySD.initsuperdropsbinary_src import read_initsuperdrops as rsupers
from pySD.gbxboundariesbinary_src import read_gbxboundaries as rgrid
from pySD.gbxboundariesbinary_src import create_gbxboundaries as cgrid

### ---------------------------------------------------------------- ###
### ----------------------- INPUT PARAMETERS ----------------------- ###
### ---------------------------------------------------------------- ###
### --- essential paths and filenames --- ###
# path and filenames for creating initial SD conditions
constsfile = path2CLEO+"/libs/cleoconstants.hpp"
binpath = path2build+"/bin/"
sharepath = path2build+"/share/"
initSDsfile = sharepath+"shima2009_dimlessSDsinit.dat"
gridfile = sharepath+"shima2009_dimlessGBxboundaries.dat"

# path and file names for plotting results
setupfile = binpath+"shima2009_setup.txt"
dataset = binpath+"shima2009_sol.zarr"

# booleans for [making, saving] initialisation figures
isfigures = [True, True]
savefigpath = path2build+"/bin/" # directory for saving figures

### --- settings for 0-D Model gridbox boundaries --- ###
zgrid = np.asarray([0, 100])
xgrid = np.asarray([0, 100])
ygrid = np.asarray([0, 100])

### --- settings for initial superdroplets --- ###
# settings for superdroplet coordinates
nsupers = {0: 4096}
coord_params = ["false"]

# settings for distirbution from exponential in droplet volume
# peak of volume exponential distribution [m]
volexpr0 = 30.531e-6
numconc = 2**(23)                     # total no. conc of real droplets [m^-3]
rspan = [1e-8, 9e-5]                # max and min range of radii to sample [m]

samplevol = rgrid.calc_domainvol(zgrid, xgrid, ygrid)
xiprobdist = probdists.VolExponential(volexpr0, rspan)
radiigen = rgens.SampleLog10RadiiGen(rspan)  # radii are sampled from rspan [m]
dryradiigen = rgens.MonoAttrGen(1e-16)       # all SDs have negligible solute [m]

coord3gen = None                        # do not generate superdroplet coords
coord1gen = None
coord2gen = None

### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###

### ---------------------------------------------------------------- ###
### ------------------- BINARY FILES GENERATION--------------------- ###
### ---------------------------------------------------------------- ###
### --- ensure build, share and bin directories exist --- ###
if path2CLEO == path2build:
  raise ValueError("build directory cannot be CLEO")
else:
  Path(path2build).mkdir(exist_ok=True)
  Path(sharepath).mkdir(exist_ok=True)
  Path(binpath).mkdir(exist_ok=True)
  if isfigures[1]:
    Path(savefigpath).mkdir(exist_ok=True)

### --- delete any existing initial conditions --- ###
os.system("rm "+gridfile)
os.system("rm "+initSDsfile)

### ----- write gridbox boundaries binary ----- ###
cgrid.write_gridboxboundaries_binary(gridfile, zgrid, xgrid,
                                     ygrid, constsfile)
rgrid.print_domain_info(constsfile, gridfile)

### ----- write initial superdroplets binary ----- ###
initattrsgen = attrsgen.AttrsGenerator(radiigen, dryradiigen, xiprobdist,
                                               coord3gen, coord1gen, coord2gen)
csupers.write_initsuperdrops_binary(initSDsfile, initattrsgen,
                                                  configfile, constsfile,
                                                  gridfile, nsupers, numconc)
rsupers.print_initSDs_infos(initSDsfile, configfile, constsfile, gridfile)

### ----- show (and save) plots of binary file data ----- ###
if isfigures[0]:
    rgrid.plot_gridboxboundaries(constsfile, gridfile,
                                 savefigpath, isfigures[1])
    rsupers.plot_initGBxs_distribs(configfile, constsfile, initSDsfile,
                                              gridfile, savefigpath, isfigures[1], "all")
    plt.close()
### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###

def run_exectuable(path2build, dataset, executable, configfile):
  ''' delete existing dataset, the run exectuable with given config file'''
  os.chdir(path2build)
  os.system('pwd')
  os.system('rm -rf '+dataset) # delete any existing dataset
  print('Executable: '+executable)
  print('Config file: '+configfile)
  os.system(executable + ' ' + configfile)

if "golovin" in kernels:
    ### ------------------------------------------------------------ ###
    ### -------------------- RUN CLEO EXECUTABLE ------------------- ###
    ### ------------------------------------------------------------ ###
    executable = path2build+'/examples/boxmodelcollisions/golovin/src/golcolls'
    run_exectuable(path2build, dataset, executable, configfile)
    ### ------------------------------------------------------------ ###
    ### ------------------------------------------------------------ ###

    ### ------------------------------------------------------------ ###
    ### ----------------------- PLOT RESULTS ----------------------- ###
    ### ------------------------------------------------------------ ###
    # read in constants and intial setup from setup .txt file
    config = pysetuptxt.get_config(setupfile, nattrs=3, isprint=True)
    consts = pysetuptxt.get_consts(setupfile, isprint=True)
    gbxs = pygbxsdat.get_gridboxes(gridfile, consts["COORD0"], isprint=True)

    time = pyzarr.get_time(dataset).secs
    sddata = pyzarr.get_supers(dataset, consts)

    # 4. plot results
    tplt = [0, 1200, 2400, 3600]
    # 0.2 factor for guassian smoothing
    smoothsig = 0.62*(config["maxnsupers"]**(-1/5))
    plotwitherr = True

    savename = savefigpath + "golovin_validation.png"
    fig, ax = shima2009fig.plot_validation_figure(plotwitherr, time,
                                        sddata, tplt, gbxs["domainvol"],
                                        numconc, volexpr0, smoothsig,
                                        savename=savename, withgol=True)
    ### ------------------------------------------------------------ ###
    ### ------------------------------------------------------------ ###

if "long" in kernels:
    ### ------------------------------------------------------------ ###
    ### -------------------- RUN CLEO EXECUTABLE ------------------- ###
    ### ------------------------------------------------------------ ###
    executable = path2build+'/examples/boxmodelcollisions/long/src/longcolls'
    run_exectuable(path2build, dataset, executable, configfile)
    ### ------------------------------------------------------------ ###
    ### ------------------------------------------------------------ ###

    ### ------------------------------------------------------------ ###
    ### ----------------------- PLOT RESULTS ----------------------- ###
    ### ------------------------------------------------------------ ###
    # read in constants and intial setup from setup .txt file
    config = pysetuptxt.get_config(setupfile, nattrs=3, isprint=True)
    consts = pysetuptxt.get_consts(setupfile, isprint=True)
    gbxs = pygbxsdat.get_gridboxes(gridfile, consts["COORD0"], isprint=True)

    time = pyzarr.get_time(dataset).secs
    sddata = pyzarr.get_supers(dataset, consts)

    # 4. plot results
    tplt = [0, 600, 1200, 1800]
    # 0.2 factor for guassian smoothing
    smoothsig = 0.62*(config["maxnsupers"]**(-1/5))
    plotwitherr = False

    savename = savefigpath + "long_validation.png"
    fig, ax = shima2009fig.plot_validation_figure(plotwitherr, time,
                                                  sddata, tplt, gbxs["domainvol"],
                                                  numconc, volexpr0, smoothsig,
                                                  savename=savename)
    ### ------------------------------------------------------------ ###
    ### ------------------------------------------------------------ ###

if "lowlist" in kernels:
    ### ------------------------------------------------------------ ###
    ### -------------------- RUN CLEO EXECUTABLE ------------------- ###
    ### ------------------------------------------------------------ ###
    executable = path2build+'/examples/boxmodelcollisions/lowlist/src/lowlistcolls'
    run_exectuable(path2build, dataset, executable, configfile)
    ### ------------------------------------------------------------ ###
    ### ------------------------------------------------------------ ###

    ### ------------------------------------------------------------ ###
    ### ----------------------- PLOT RESULTS ----------------------- ###
    ### ------------------------------------------------------------ ###
    # read in constants and intial setup from setup .txt file
    config = pysetuptxt.get_config(setupfile, nattrs=3, isprint=True)
    consts = pysetuptxt.get_consts(setupfile, isprint=True)
    gbxs = pygbxsdat.get_gridboxes(gridfile, consts["COORD0"], isprint=True)

    time = pyzarr.get_time(dataset).secs
    sddata = pyzarr.get_supers(dataset, consts)

    # 4. plot results
    tplt = [0, 600, 1200, 1800, 2400, 3600]
    # 0.2 factor for guassian smoothing
    smoothsig = 0.62*(config["maxnsupers"]**(-1/5))
    plotwitherr = False

    savename = savefigpath + "lowlist_validation.png"
    fig, ax = shima2009fig.plot_validation_figure(plotwitherr, time,
                                                  sddata, tplt, gbxs["domainvol"],
                                                  numconc, volexpr0, smoothsig,
                                                  savename=savename)
    ### ------------------------------------------------------------ ###
    ### ------------------------------------------------------------ ###
