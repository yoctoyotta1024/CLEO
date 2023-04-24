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
# gen = thermogen.ConstUniformThermo(P_INIT, TEMP_INIT, None,
#                                        qc_init, W_INIT, U_INIT, V_INIT,
#                                        relh=relh_init, constsfile)

PRESS0 = 101500 # [Pa]
THETA = 289 # [K]
qvap = 0.0075 # [Kg/Kg]
qcond = 0.0 # [Kg/Kg]
WMAX = 0.6 # [m/s]
VVEL = None # [m/s]
Zlength = 1500 # [m]
Xlength = 1500 # [m]
inputs = cthermo.thermoinputsdict(configfile, constsfile)

# gen = thermogen.ConstHydrostaticAdiabat(PRESS0, THETA, qvap, qcond, WMAX, 
#                                         Zlength, Xlength, VVEL,
#                                         inputs["G"], inputs["CP_DRY"],
#                                         inputs["RGAS_DRY"], inputs["RGAS_V"])
gen = thermogen.ConstThermo2Dflowfield(PRESS0, THETA, "sratio", qcond, WMAX,
                                        Zlength, Xlength, VVEL,
                                        qparam=1.05, constsfile=constsfile)
cthermo.write_thermodynamics_binary(thermofile, gen, configfile,
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

        if "wvel" in thermodata.keys():
            self.wvel = np.reshape(thermodata["wvel"], reshape) 
            if "uvel" in thermodata.keys():
                self.uvel = np.reshape(thermodata["uvel"], reshape) 
                if "vvel" in thermodata.keys():
                    self.vvel = np.reshape(thermodata["vvel"], reshape)  

    def meanytime(self, var):

        return np.mean(var, axis=(0,1))

    def meanxytime(self, var):

        return np.mean(var, axis=(0,1, 2))

ngridboxes=int(np.prod(ndims))
thermodata = rthermo.get_thermodynamics_from_thermofile(thermofile, ngridboxes,
                                               inputs=inputs)
d = Thermo2Davg(thermodata, ndims, inputs["ntime"])

meanpress = d.meanytime(d.press)
meantemp = d.meanytime(d.temp)
meanqvap = d.meanytime(d.qvap)
meanqcond = d.meanytime(d.qcond)
meanwvel = d.meanytime(d.wvel)
meanuvel = d.meanytime(d.uvel)

xxh, zzh = np.meshgrid(xhalf, zhalf, indexing="ij") # dims [xdims, zdims]
xxf, zzf = np.meshgrid(xfull, zfull, indexing="ij") # dims [xdims, zdims]

meanrelh, meansupersat = rthermo.relative_humidity(meanpress, meantemp,
                                                meanqvap, inputs["Mr_ratio"])

fig, axs = plt.subplots(nrows=2, ncols=2)
axs = axs.flatten()
axs[0].pcolormesh(xxh[:,:], zzh[:,:], meanpress, cmap="plasma_r")
axs[1].pcolormesh(xxh[:,:], zzh[:,:], meantemp, cmap="plasma_r")
axs[2].pcolormesh(xxh[:,:], zzh[:,:], meanqvap, cmap="plasma_r")
#axs[3].pcolormesh(xxh[:,:], zzh[:,:], meanqcond, cmap="plasma_r")
axs[3].pcolormesh(xxh[:,:], zzh[:,:], meansupersat, cmap="plasma_r")
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

fig, axs = plt.subplots(nrows=2, ncols=2)
axs = axs.flatten()
labs = ["P /mbar", "T /K", "q$_{vapour}$", "supersaturation"]
legs = []
axs[0].plot(d.meanxytime(d.press), zfull)
legs.append("<P>={:.1f}Pa".format(np.mean(d.meanxytime(d.press))))
axs[1].plot(d.meanxytime(d.temp), zfull)
legs.append("<T>={:.1f}K".format(np.mean(d.meanxytime(d.temp))))
axs[2].plot(d.meanxytime(d.qvap), zfull)
legs.append("<q$_{v}$>="+"{:.3g}".format(np.mean(d.meanxytime(d.qvap))))
axs[3].plot(np.mean(meansupersat, axis=0), zfull)
legs.append("<s>="+"{:.3g}".format(np.mean(meansupersat)))

axs[3].vlines(0.0, np.amin(zfull), np.amax(zfull),
              color="grey", linestyle="--")
s0 = zfull[np.argmin(abs(np.mean(meansupersat, axis=0)))]
stext = "s$_{0}$="+"{:.2f}m".format(s0)
axs[3].text(0.95, 0.05, stext, transform=axs[3].transAxes,
            horizontalalignment="right")

for a, ax in enumerate(axs):
  ax.set_xlabel(labs[a])
  ax.set_ylabel("z / m")
  ax.text(0.95, 0.85, legs[a], transform=ax.transAxes, horizontalalignment="right")
fig.tight_layout()
if isfigures[1]:
    fig.savefig(binpath+"thermoprofile1D.png", dpi=400,
                bbox_inches="tight", facecolor='w', format="png")
    print("Figure .png saved as: "+binpath+"thermoprofile1D.png")
plt.show()