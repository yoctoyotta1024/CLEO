'''
----- CLEO -----
File: shima2009.py
Project: golovin
Created Date: Friday 17th November 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
Last Modified: Friday 17th November 2023
Modified By: CB
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
Copyright (c) 2023 MPI-M, Clara Bayley
-----
File Description:
Script compiles and runs CLEO adia0D to create the
data and plots as in Shima et al. 2009 to show
comparision of SDM 0D box model of collision-
coalescence with Golovin's analytical solution
'''

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

path2CLEO = sys.argv[1]
path2build = sys.argv[2]
configfile = sys.argv[3]

sys.path.append(path2CLEO)  # for imports from pySD package
sys.path.append(path2CLEO+"/examples/exampleplotting/") # for imports from example plotting package

from plotssrc import shima2009fig
from pySD.output_src import *
from pySD.initsuperdropsbinary_src import *
from pySD.gbxboundariesbinary_src import read_gbxboundaries as rgrid
from pySD.gbxboundariesbinary_src import create_gbxboundaries as cgrid

############### INPUTS ##################
# path and filenames for creating SD initial conditions and for running model
constsfile = path2CLEO+"/libs/cleoconstants.hpp"
binpath = path2build+"/bin/"
sharepath = path2build+"/share/"
initSDsfile = sharepath+"golovin_dimlessSDsinit.dat"
gridfile = sharepath+"golovin_dimlessGBxboundaries.dat"

# path and file names for plotting results
setupfile = binpath+"golovin_setup.txt"
dataset = binpath+"golovin_sol.zarr"

# booleans for [making, showing] initialisation figures
isfigures = [True, False]

# settings for 0D Model (number of SD and grid coordinates)
nsupers = {0: 2048}
coord_params = ["false"]
zgrid = np.asarray([0, 100])
xgrid = np.asarray([0, 100])
ygrid = np.asarray([0, 100])

# settings for distirbution from exponential in droplet volume
# peak of volume exponential distribution [m]
volexpr0 = 30.531e-6
numconc = 2**(23)                     # total no. conc of real droplets [m^-3]
rspan = [1e-8, 9e-5]                # max and min range of radii to sample [m]
randomr = True                        # sample radii range randomly or not

samplevol = rgrid.calc_domainvol(zgrid, xgrid, ygrid)
radiiprobdist = radiiprobdistribs.VolExponential(volexpr0, rspan)
radiigen = initattributes.SampleDryradiiGen(
    rspan, randomr)  # radii are sampled from rspan [m]
coord3gen = None                        # do not generate superdroplet coords
coord1gen = None
coord2gen = None

# ### 1. create files for gridbox boundaries and initial SD conditions
# os.system("rm "+gridfile)
# os.system("rm "+initSDsfile)
# cgrid.write_gridboxboundaries_binary(gridfile, zgrid, xgrid,
#                                      ygrid, constsfile)
# rgrid.print_domain_info(constsfile, gridfile)

# initattrsgen = initattributes.InitManyAttrsGen(radiigen, radiiprobdist,
#                                                coord3gen, coord1gen, coord2gen)
# create_initsuperdrops.write_initsuperdrops_binary(initSDsfile, initattrsgen,
#                                                   configfile, constsfile,
#                                                   gridfile, nsupers, numconc)
# read_initsuperdrops.print_initSDs_infos(initSDsfile, configfile, constsfile, gridfile)

# if isfigures[0]:
#     rgrid.plot_gridboxboundaries(constsfile, gridfile,
#                                  binpath, isfigures[1])
#     read_initsuperdrops.plot_initGBxsdistribs(configfile, constsfile, initSDsfile,
#                                               gridfile, binpath, isfigures[1], "all")
# plt.close()

# # 2. compile and the run model
# os.chdir(path2build)
# os.system('pwd')
# os.system('rm -rf '+dataset)
# os.system("make clean && make -j 64 gol0D")
# executable = path2build+"/examples/golovin/src/gol0D"
# os.system(executable + ' ' + configfile)

# 3. load and plot results
# read in constants and intial setup from setup .txt file
config = pysetuptxt.get_config(setupfile, nattrs=3, isprint=True)
consts = pysetuptxt.get_consts(setupfile, isprint=True)
gbxs = pygbxsdat.get_gridboxes(gridfile, consts["COORD0"], isprint=True)

time = pyzarr.get_time(dataset).secs
sddata = pyzarr.get_supers(dataset, consts)

# 4. plot results
tplt = [0, 1200, 2400, 3600]
# 0.2 factor for guassian smoothing
smoothsig = 0.62*(config["totnsupers"]**(-1/5))
plotwitherr = True

savename = binpath + "golovin_validation.png"
fig, ax = shima2009fig.golovin_validation_figure(plotwitherr, time,
                                    sddata, tplt, gbxs["domainvol"],
                                    numconc, volexpr0, smoothsig,
                                    savename=savename)