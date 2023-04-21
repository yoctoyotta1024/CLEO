import numpy as np
import matplotlib.pyplot as plt 

import pySD.gbxboundariesbinary_src.read_gbxboundaries as rgrid

abspath = "/Users/yoctoyotta1024/Documents/b1_springsummer2023/CLEO/"
constsfile = abspath+"libs/claras_SDconstants.hpp"
configfile = abspath+"src/config/config.txt"

gridfile = abspath+"build/share/dimlessGBxboundaries.dat"
binpath = abspath+"build/bin/"

COORD0 = rgrid.get_COORD0_from_constsfile(constsfile)
gbxbounds, ndims = rgrid.read_dimless_gbxboundaries_binary(gridfile,
                                                            COORD0=COORD0,
                                                            return_ndims=True)
zhalf, xhalf, yhalf = rgrid.halfcoords_from_gbxbounds(gbxbounds)
zfull, xfull, yfull = rgrid.fullcell_fromhalfcoords(zhalf, xhalf, yhalf)
print(ndims, "ie. shapes should be: ", np.flip(ndims))
print(zfull)
print(xfull)
print(yfull)
print("----")

yy, xx, zz = np.meshgrid(yfull, xfull, zfull, indexing="ij") # dims [xdims, zdims]
print(xx.shape)

zfulls, xfulls, yfulls = rgrid.fullcoords_forallgridboxes(gbxbounds, ndims)

for j in range(ndims[2]):
    for i in range(ndims[1]):
        for k in range(ndims[0]):
            ii = int(j*ndims[0]*ndims[1] + i*ndims[0] + k)
            zxyfull = [zfulls[ii], xfulls[ii], yfulls[ii]]
            zxymesh = [zz[j,i,k], xx[j,i,k], yy[j,i,k]]
            print("GBx: ", ii, (k,i,j), "-> coord centre: ", zxyfull,
                   " ie. ", zxymesh, zxyfull==zxymesh)
                  
print("zs", np.all(zfulls == zz.flatten()))
print("xs", np.all(xfulls == xx.flatten()))
print("ys", np.all(yfulls == yy.flatten()))

shape = np.flip(ndims)
zmesh = np.reshape(zfulls, shape)
xmesh = np.reshape(xfulls, shape)
ymesh = np.reshape(yfulls, shape)
print(np.all(zmesh == zz))
print(np.all(xmesh == xx))
print(np.all(ymesh == yy))
