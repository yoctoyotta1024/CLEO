'''
----- CLEO -----
File: speedtest.py
Project: speedtest
Created Date: Friday 17th November 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
Last Modified: Saturday 9th December 2023
Modified By: CB
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
Copyright (c) 2023 MPI-M, Clara Bayley
-----
File Description:
Script compiles and runs CLEO speedtest to 
check performance of CLEO usign different
build configurations (e.g. serial, OpenmP
and CUDA parallelism).
'''

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

path2CLEO = sys.argv[1]
path2build = sys.argv[2]
configfile = sys.argv[3]
outputdir = sys.argv[4]
buildtype = sys.argv[5]

sys.path.append(path2CLEO)  # for imports from pySD package
sys.path.append(path2CLEO+"/examples/exampleplotting/") # for imports from example plotting package

from plotssrc import pltsds, pltmoms
from pySD.sdmout_src import *
from pySD.initsuperdropsbinary_src import *
from pySD.gbxboundariesbinary_src import read_gbxboundaries as rgrid
from pySD.gbxboundariesbinary_src import create_gbxboundaries as cgrid
from pySD.initsuperdropsbinary_src import initattributes as iattrs
from pySD.initsuperdropsbinary_src import radiiprobdistribs as rprobs 
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
constsfile = path2CLEO+"/libs/cleoconstants.hpp"
binpath = path2build+"/bin/"
sharepath = path2build+"/share/"
gridfile = sharepath+"spd_dimlessGBxboundaries.dat"
initSDsfile = sharepath+"spd_dimlessSDsinit.dat"
thermofile =  sharepath+"/spd_dimlessthermo.dat"

# path and file names for plotting results
setupfile = binpath+"spd_setup.txt"
statsfile = binpath+"spd_stats.txt"
dataset = binpath+"spd_sol.zarr"

### --- plotting initialisation figures --- ###
isfigures = [True, True] # booleans for [making, saving] initialisation figures
savefigpath = outputdir   # directory for saving figures
SDgbxs2plt = [0] # gbxindex of SDs to plot (nb. "all" can be very slow)
outdatafile = outputdir+"spd_allstats.txt" # file to write out stats to

### --- settings for 2-D gridbox boundaries --- ###
zgrid = [0, 1500, 50]     # evenly spaced zhalf coords [zmin, zmax, zdelta] [m]
xgrid = [0, 1500, 50]     # evenly spaced xhalf coords [m]
ygrid = np.array([0, 25, 50])  # array of yhalf coords [m]

### --- settings for initial superdroplets --- ###
# settings for initial superdroplet coordinates
zlim = 1500       # max z coord of superdroplets
npergbx = 8       # number of superdroplets per gridbox 

# [min, max] range of initial superdroplet radii (and implicitly solute masses)
rspan                = [3e-9, 3e-6] # [m]

# settings for initial superdroplet multiplicies
# (from bimodal Lognormal distribution)
geomeans             = [0.02e-6, 0.15e-6]               
geosigs              = [1.4, 1.6]                    
scalefacs            = [6e6, 4e6]   
numconc = np.sum(scalefacs)

### --- settings for 3D Thermodyanmics --- ###
PRESS0 = 100000 # [Pa]
THETA = 298.15  # [K]
qcond = 0.0     # [Kg/Kg]
WMAX = 3.0      # [m/s]
VVEL = 1.0      # [m/s]
Zlength = 1500  # [m]
Xlength = 1500  # [m]
qvapmethod = "sratio"
Zbase = 750     # [m]
sratios = [0.85, 1.1] # s_ratio [below, above] Zbase
### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###

### ---------------------------------------------------------------- ###
### --------------------- FUNCTION DEFINITIONS --------------------- ###
### ---------------------------------------------------------------- ###
def read_statsfile(statsfile):

  stats = {}
  with open(statsfile, 'r') as file:
    for line in file:
            # Check if the line starts with '###'
            if not line.startswith('###'):
                # Process the line
                line = line.strip().split()
                stats[line[0]] = float(line[1])
  
  return stats

def write_outstats(outdatafile, stats):
  ''' if outdatafile doesn't already exist, creates new file with 
  a header. else appends to end of file '''
  
  try:
    # Try to open the file for exclusive creation
    with open(outdatafile, 'x') as file:
        # Perform operations on the new file if needed
        header = "### Wall Clock time For Timestepping\n"
        header += "### columns are: "
        header += "Test Number gpus_cpus/s cpus/s serial/s\n"
        file.write(header)
    print(f"--- new stats output: '{outdatafile}' created ---")
  except FileExistsError:
    print(f"stats output file '{outdatafile}' already exists")

  # count existing number of lines in 
  with open(outdatafile, 'r') as file:
    lines = file.readlines()
    headerlines = [line for line in lines if line.startswith("###")]
    nlines = len(lines)
    nheader = len(headerlines)
        
  # write stats["tstep"] value to output file
  testnumber = int(nlines - nheader)
  with open(outdatafile, 'a') as file:
    line = "\n"+str(testnumber)+" "+str(stats["tstep"])
    file.write(line)
   
### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###

# ### ---------------------------------------------------------------- ###
# ### ------------------- BINARY FILES GENERATION--------------------- ###
# ### ---------------------------------------------------------------- ###
# ### --- ensure build, share and bin directories exist --- ###
# if path2CLEO == path2build:
#   raise ValueError("build directory cannot be CLEO")
# else:
#   Path(path2build).mkdir(exist_ok=True) 
#   Path(sharepath).mkdir(exist_ok=True) 
#   Path(binpath).mkdir(exist_ok=True) 
#   if isfigures[1]:
#     Path(savefigpath).mkdir(exist_ok=True) 
# os.system("rm "+gridfile)
# os.system("rm "+initSDsfile)
# os.system("rm "+thermofile[:-4]+"*")

# ### ----- write gridbox boundaries binary ----- ###
# cgrid.write_gridboxboundaries_binary(gridfile, zgrid, xgrid, ygrid, constsfile)
# rgrid.print_domain_info(constsfile, gridfile)

# ### ----- write thermodyanmics binaries ----- ###
# thermodyngen = thermogen.SimpleThermo2Dflowfield(configfile, constsfile, PRESS0,
#                                                 THETA, qvapmethod, sratios, Zbase,
#                                                 qcond, WMAX, Zlength, Xlength,
#                                                 VVEL)
# cthermo.write_thermodynamics_binary(thermofile, thermodyngen, configfile,
#                                     constsfile, gridfile)


# ### ----- write initial superdroplets binary ----- ###
# nsupers = iattrs.nsupers_at_domain_base(gridfile, constsfile, npergbx, zlim)
# coord3gen = iattrs.SampleCoordGen(True) # sample coord3 randomly
# coord1gen = iattrs.SampleCoordGen(True) # sample coord1 randomly
# coord2gen = iattrs.SampleCoordGen(True) # sample coord2 randomly 
# radiiprobdist = rprobs.LnNormal(geomeans, geosigs, scalefacs)
# radiigen = iattrs.SampleDryradiiGen(rspan, True) # randomly sample radii from rspan [m]

# initattrsgen = iattrs.InitManyAttrsGen(radiigen, radiiprobdist,
#                                         coord3gen, coord1gen, coord2gen)
# csupers.write_initsuperdrops_binary(initSDsfile, initattrsgen, 
#                                       configfile, constsfile,
#                                       gridfile, nsupers, numconc)

# ### ----- show (and save) plots of binary file data ----- ###
# if isfigures[0]:
#   rgrid.plot_gridboxboundaries(constsfile, gridfile,
#                                savefigpath, isfigures[1])
#   rthermo.plot_thermodynamics(constsfile, configfile, gridfile,
#                               thermofile, savefigpath, isfigures[1])
#   rsupers.plot_initGBxsdistribs(configfile, constsfile, initSDsfile,
#                               gridfile, savefigpath, isfigures[1],
#                               SDgbxs2plt) 
#   plt.close()
# ### ---------------------------------------------------------------- ###
# ### ---------------------------------------------------------------- ###

# ### ---------------------------------------------------------------- ###
# ### -------------------- COMPILE AND RUN CLEO ---------------------- ###
# ### ---------------------------------------------------------------- ###
# os.chdir(path2build)
# os.system('pwd')
# os.system('rm -rf '+dataset)
# os.system('make clean && make -j 64 spdtest')
# executable = path2build+'/examples/speedtest/src/spdtest'
# os.system(executable + ' ' + configfile)
# ### ---------------------------------------------------------------- ###
# ### ---------------------------------------------------------------- ###

### ---------------------------------------------------------------- ###
### ------------------- WRITE RESULTS TO FILE ---------------------- ###
### ---------------------------------------------------------------- ###
print("--- reading runtime statistics ---")     
stats = read_statsfile(statsfile)
for key, value in stats.items():
    print(key+": {:.3f}s".format(value))
write_outstats(outdatafile, stats)
print("--- runtime stats written to file ---")
### ---------------------------------------------------------------- ###
### ---------------------------------------------------------------- ###
