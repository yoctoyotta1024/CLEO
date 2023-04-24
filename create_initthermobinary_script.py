import numpy as np
from pySD.thermobinary_src import thermogen
from pySD.thermobinary_src import create_thermodynamics as cthermo
from pySD.thermobinary_src import read_thermodynamics as rthermo

# abspath = "/Users/yoctoyotta1024/Documents/b1_springsummer2023/CLEO/"
abspath = "/home/m/m300950/CLEO/"
#abspath = sys.argv[1]
constsfile = abspath+"libs/claras_SDconstants.hpp"
configfile = abspath+"src/config/config.txt"

spath = abspath+"build/share/"
gridfile = spath+"dimlessGBxboundaries.dat"
thermofile =  spath+"dimlessthermodynamics.dat"
binpath = abspath+"build/bin/"

### booleans for [making+showing, saving] figures
isfigures = [True, True]

### initial thermodynamic conditions for all gridboxes ###
# P_INIT = 100000.0                       # initial pressure [Pa]
# TEMP_INIT = 273.15                      # initial parcel temperature [T]
# relh_init = 95.0                        # initial relative humidity (%)
# qc_init = 0.0                           # initial liquid water content []
# W_INIT = 0.0                            # initial vertical (z) velocity [m/s]
# U_INIT = 0.0                           # initial horizontal x velocity [m/s]
# V_INIT = 0.0                             # initial horizontal y velocity [m/s]
# gen = thermogen.ConstUniformThermo(P_INIT, TEMP_INIT, relh_init,
#                                        qc_init, W_INIT, U_INIT, V_INIT,
#                                        constsfile)

PRESS0 = 101500 # [Pa]
THETA = 289 # [K]
qvap = 0.0075 # [Kg/Kg]
qcond = 0.0 # [Kg/Kg]
WMAX = 0.6 # [m/s]
VVEL = None # [m/s]
Zlength = 1500 # [m]
Xlength = 1500 # [m]
inputs = cthermo.thermoinputsdict(configfile, constsfile)
gen = thermogen.ConstHydrostaticAdiabat(PRESS0, THETA, qvap, qcond, WMAX, 
                                        Zlength, Xlength, VVEL,
                                        inputs["G"], inputs["CP_DRY"],
                                        inputs["RGAS_DRY"], inputs["RGAS_V"])

thermodata = cthermo.write_thermodynamics_binary(thermofile, gen, configfile,
                                    constsfile, gridfile)

# if isfigures[0]:
#     #t2plt = 0.0
#     t2plt =  "all"
#     rthermo.plot_thermodynamics_timeslice(constsfile, configfile, gridfile,
#                                           thermofile, t2plt, binpath,
#                                           isfigures[1])

from pySD.gbxboundariesbinary_src import read_gbxboundaries as rgrid 
import matplotlib.pyplot as plt

gbxbounds, ndims = rgrid.read_dimless_gbxboundaries_binary(gridfile,
                                                            COORD0=inputs["COORD0"],
                                                            return_ndims=True)
zhalf, xhalf, yhalf = rgrid.halfcoords_from_gbxbounds(gbxbounds)
zfull, xfull, yfull = rgrid.fullcell_fromhalfcoords(zhalf, xhalf, yhalf)


class Thermo2Davg:

    def __init__(self, thermodata, ndims, ntime):
        
        reshape = [ntime] + list(np.flip(ndims))
        
        self.press = np.reshape(thermodata["press"], reshape)
        self.temp = np.reshape(thermodata["temp"], reshape)
        self.qvap= np.reshape(thermodata["qvap"], reshape)
        self.qcond = np.reshape(thermodata["qcond"], reshape)
        self.wvel = np.array([])
        self.uvel = np.array([])
        self.vvel = np.array([])

        if any(thermodata["wvel"]):
            self.wvel = np.reshape(thermodata["wvel"], reshape) 
            if any(thermodata["uvel"]):
                self.uvel = np.reshape(thermodata["uvel"], reshape) 
                if any(thermodata["vvel"]):
                    self.vvel = np.reshape(thermodata["vvel"], reshape)  

    def meanYtime(self, var):

        return np.mean(var, axis=(0,1))

# redim = cthermo.DimlessThermodynamics(inputs=inputs)
# thermodata = redim.redimensionalise(thermodata)
d = Thermo2Davg(thermodata, ndims, inputs["ntime"])

meanpress = d.meanYtime(d.press)
meantemp = d.meanYtime(d.temp)
meanqvap = d.meanYtime(d.qvap)
meanqcond = d.meanYtime(d.qcond)
meanwvel = d.meanYtime(d.wvel)
meanuvel = d.meanYtime(d.uvel)

xxh, zzh = np.meshgrid(xhalf, zhalf, indexing="ij") # dims [xdims, zdims]
xxf, zzf = np.meshgrid(xfull, zfull, indexing="ij") # dims [xdims, zdims]

fig, axs = plt.subplots(nrows=2, ncols=2)
axs = axs.flatten()
axs[0].pcolormesh(xxh[:,:], zzh[:,:], meanpress, cmap="plasma_r")
axs[1].pcolormesh(xxh[:,:], zzh[:,:], meantemp, cmap="plasma_r")
axs[2].pcolormesh(xxh[:,:], zzh[:,:], meanqvap, cmap="plasma_r")
axs[3].pcolormesh(xxh[:,:], zzh[:,:], meanqcond, cmap="plasma_r")
fig.tight_layout()
plt.show()

fig, axs = plt.subplots(3)
for i in range(ndims[1]):
    axs[0].plot(meanwvel[i,:], zfull)
    axs[0].set_ylabel("z")
for j in range(ndims[0]):
    axs[1].plot(xfull, meanwvel[:,i])
    axs[1].set_xlabel("x")
axs[2].pcolormesh(xxh[:,:], zzh[:,:], meanwvel, cmap="plasma_r")
fig.tight_layout()
plt.show()

fig, axs = plt.subplots(3)
for i in range(ndims[1]):
    axs[0].plot(meanuvel[i,:], zfull)
    axs[0].set_ylabel("z")
for j in range(ndims[0]):
    axs[1].plot(xfull, meanuvel[:,i])
    axs[1].set_xlabel("x")
axs[2].pcolormesh(xxh[:,:], zzh[:,:], meanuvel, cmap="plasma_r")
fig.tight_layout()
plt.show()

fig, ax = plt.subplots()
print(xxf.shape, zzf.shape, meanuvel.shape, meanwvel.shape)
ax.quiver(xxf, zzf, meanuvel, meanwvel)
plt.show()