import numpy as np
from pySD.thermobinary_src import thermogen
from pySD.thermobinary_src import create_thermodynamics as cthermo
from pySD.thermobinary_src import read_thermodynamics as rthermo

abspath = "/Users/yoctoyotta1024/Documents/b1_springsummer2023/CLEO/"
constsfile = abspath+"libs/claras_SDconstants.hpp"
configfile = abspath+"src/config/config.txt"

spath = abspath+"build/share/"
gridfile = spath+"dimlessGBxboundaries.dat"
thermofile =  spath+"dimlessthermodynamics.dat"
binpath = abspath+"build/bin/"

### booleans for [making+showing, saving] figures
isfigures = [True, True]

### initial thermodynamic conditions for all gridboxes ###
P_INIT = 100000.0                       # initial pressure [Pa]
TEMP_INIT = 273.15                      # initial parcel temperature [T]
relh_init = 95.0                        # initial relative humidity (%)
qc_init = 0.0                           # initial liquid water content []
W_INIT = 1.0                            # initial vertical (z) velocity [m/s]
U_INIT = 0.5                           # initial horizontal x velocity [m/s]
V_INIT = 3.0                             # initial horizontal y velocity [m/s]

# gen = thermogen.ConstUniformThermo(P_INIT, TEMP_INIT, relh_init,
#                                        qc_init, W_INIT, U_INIT, V_INIT,
#                                        constsfile)
# cthermo.write_thermodynamics_binary(thermofile, gen, configfile,
#                                     constsfile, gridfile)

# if isfigures[0]:
#     #t2plt = 0.0
#     t2plt =  "all"
#     rthermo.plot_thermodynamics_timeslice(constsfile, configfile, gridfile,
#                                           thermofile, t2plt, binpath,
#                                           isfigures[1])



import pySD.gbxboundariesbinary_src.read_gbxboundaries as rgrid
import matplotlib.pyplot as plt 
inputs = cthermo.thermoinputsdict(configfile, constsfile)
gbxbounds, ndims = rgrid.read_dimless_gbxboundaries_binary(gridfile,
                                                COORD0=inputs["COORD0"],
                                                return_ndims=True)
zhalf, xhalf, yhalf = rgrid.halfcoords_from_gbxbounds(gbxbounds)
zfull, xfull, yfull = rgrid.get_fullcoords_from_gridfile(gridfile, COORD0=inputs["COORD0"])
shape = list(np.flip(ndims))
ntime = 81
print(ndims, np.flip(ndims))
print(zfull)
print(xfull)
print(yfull)
print("----")

yy, xx, zz = np.meshgrid(yfull, xfull, zfull, indexing="ij") # dims [xdims, zdims]
print(xx.shape)

zfulls = np.tile(zfull, int(ndims[1]*ndims[2])) # zfull of every gridbox in order of gbxindex
xfulls = np.tile(np.repeat(xfull, ndims[0]), int(ndims[2]))
yfulls = np.repeat(yfull, ndims[0]*ndims[1])

for j in range(ndims[2]):
    for i in range(ndims[1]):
        for k in range(ndims[0]):
            ii = int(j*ndims[0]*ndims[1] + i*ndims[0] + k)
            zxyfull = [zfulls[ii], xfulls[ii], yfulls[ii]]
            zxymesh = [zz[j,i,k], xx[j,i,k], yy[j,i,k]]
            print("GBx: ", ii, [k,i,j], "->", zxyfull,
                   " ie. ", zxymesh, zxyfull==zxymesh)
                  
print("zs", np.all(zfulls == zz.flatten()))
print("xs", np.all(xfulls == xx.flatten()))
print("ys", np.all(yfulls == yy.flatten()))

zmesh = np.reshape(zfulls, shape)
xmesh = np.reshape(xfulls, shape)
ymesh = np.reshape(yfulls, shape)
print(np.all(zmesh == zz))
print(np.all(xmesh == xx))
print(np.all(ymesh == yy))

def hydrostatic_pressure_profile(z, ntime):

    p0 = 1000 #hPa
    z0 = 1250 #km
    p = p0 * np.exp(1-z/z0)

    return np.tile(p, ntime)

press = hydrostatic_pressure_profile(zfulls, ntime)
press = np.reshape(press, [ntime]+shape)
print(press, press.shape)
meanpress = np.mean(press, axis=(0,1)) # avg over y dim and time

xxh, zzh = np.meshgrid(xhalf, zhalf, indexing="ij") # dims [xdims, zdims]

fig, ax = plt.subplots()
ax.pcolormesh(xxh[:,:], zzh[:,:], meanpress, cmap="plasma_r")

fig.tight_layout()
plt.show()

def windfield2D_wvel(x, z, ntime):

    w0 = 1
    x0 = 150/(2*np.pi)
    z0 = 1500/(np.pi)
    w = - w0 * np.sin(x/x0) * np.sin(z/z0) * 1/x0

    return np.tile(w, ntime)

def windfield2D_uvel(x, z, ntime):

    u0 = 1
    x0 = 150/(2*np.pi)
    z0 = 1500/(np.pi)
    w = - u0 * np.cos(x/x0) * np.cos(z/z0) * 1/z0

    return np.tile(w, ntime)

wvel = windfield2D_wvel(xfulls, zfulls, ntime)
wvel = np.reshape(wvel, [ntime]+shape)
wvelmean = np.mean(wvel, axis=(0,1))
fig, axs = plt.subplots(3)
for i in range(ndims[1]):
    axs[0].plot(wvelmean[i,:], zfull)
    axs[0].set_ylabel("z")
for j in range(ndims[0]):
    axs[1].plot(xfull, wvelmean[:,i])
    axs[1].set_xlabel("x")
axs[2].pcolormesh(xxh[:,:], zzh[:,:], wvelmean, cmap="plasma_r")
fig.tight_layout()
plt.show()

uvel = windfield2D_uvel(xfulls, zfulls, ntime)
uvel = np.reshape(uvel, [ntime]+shape)
uvelmean = np.mean(uvel, axis=(0,1))
fig, axs = plt.subplots(3)
for i in range(ndims[1]):
    axs[0].plot(uvelmean[i,:], zfull)
    axs[0].set_ylabel("z")
for j in range(ndims[0]):
    axs[1].plot(xfull, uvelmean[:,i])
    axs[1].set_xlabel("x")
axs[2].pcolormesh(xxh[:,:], zzh[:,:], uvelmean, cmap="plasma_r")
fig.tight_layout()
plt.show()

xxf, zzf = np.meshgrid(xfull, zfull, indexing="ij") # dims [xdims, zdims]
fig, ax = plt.subplots()
print(xxf.shape, zzf.shape, uvelmean.shape, wvelmean.shape)
ax.quiver(xxf, zzf, uvelmean, wvelmean)
plt.show()