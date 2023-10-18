'''
----- CLEO -----
File: read_thermodynamics.py
Project: thermobinary_src
Created Date: Friday 13th October 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
Last Modified: Wednesday 18th October 2023
Modified By: CB
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
Copyright (c) 2023 MPI-M, Clara Bayley
-----
File Description:
'''


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors

from .create_thermodynamics import thermoinputsdict, DimlessThermodynamics
from .thermogen import saturation_press
from ..readbinary import readbinary
from ..gbxboundariesbinary_src import read_gbxboundaries as rgrid 

def relative_humidity(press, temp, qvap, Mr_ratio):

    pv = qvap*press/(Mr_ratio + qvap) # vapour pressure
    psat = saturation_press(temp)
    relh = pv/psat
    
    qsat = Mr_ratio * psat/(press-pv) 
    supersat = qvap/qsat - 1  

    return relh, supersat

class ThermoOnGrid:

  def __init__(self, thermodata, inputs, ndims):

    dth = DimlessThermodynamics(inputs=inputs)
    thermodata = dth.redimensionalise(thermodata)
    
    ntime = inputs["ntime"]
    ndims = [int(n) for n in ndims]
    cen = [inputs["ntime"]] + list(np.flip(ndims))
    zface = [ntime, ndims[2], ndims[1], ndims[0]+1]
    xface = [ntime, ndims[2], ndims[1]+1, ndims[0]]
    yface = [ntime, ndims[2]+1, ndims[1], ndims[0]]

    self.vars = ["press", "temp", "qvap", "qcond"]
    self.press = np.reshape(thermodata["press"], cen) 
    self.temp = np.reshape(thermodata["temp"], cen) 
    self.qvap= np.reshape(thermodata["qvap"], cen) 
    self.qcond = np.reshape(thermodata["qcond"], cen) 

    if "wvel" in thermodata.keys():
        self.wvel = np.reshape(thermodata["wvel"], zface) 
        self.wvel_cens = (self.wvel[:,:,:,1:] + self.wvel[:,:,:,:-1])/2
        self.vars += ["wvel", "wvel_cens"]
        if "uvel" in thermodata.keys():
            self.uvel = np.reshape(thermodata["uvel"], xface) 
            self.uvel_cens = (self.uvel[:,:,1:,:] + self.uvel[:,:,:-1,:])/2
            self.vars += ["uvel", "uvel_cens"]
            if "vvel" in thermodata.keys():
                self.vvel = np.reshape(thermodata["vvel"], yface)  
                self.vvel_cens = (self.vvel[:,1:,:,:] + self.vvel[:,:-1,:,:])/2
                self.vars += ["vvel", "vvel_cens"]

  def __getitem__(self, key):
    if key not in self.vars:
      err = "no variable "+key+" in ThermoOnGrid"
      raise ValueError(err)         
    elif key == "press":
      return self.press
    elif key == "temp":
      return self.temp
    elif key == "qvap":
      return self.qvap
    elif key == "qcond":
      return self.qcond
    
    elif key == "wvel":
      return self.wvel
    elif key == "uvel":
      return self.uvel
    elif key == "vvel":
      return self.vvel

    elif key == "wvel_cens":
      return self.wvel_cens
    elif key == "uvel_cens":
      return self.uvel_cens
    elif key == "vvel_cens":
      return self.vvel_cens
    
  def xymean(self, var):
    '''mean over x and y'''

    return np.mean(var, axis=(1,2)) #d ims [time, z]
    
  def ytmean(self, var):
    '''mean over time and y'''

    return np.mean(var, axis=(0,1)) # dims [x, z]

  def xytmean(self, var):
    '''mean over time, x and y'''

    return np.mean(var, axis=(0,1,2)) # dims [z]

def thermovar_from_binary(var, thermofile, shape,
                          ntime, ndims, dtype,
                          isprint=True):
  
  idot = [i for i, ltr in enumerate(thermofile) if ltr == "."][-1]
  filestem, filetype = thermofile[:idot], thermofile[idot:]
  filename = filestem+"_"+var+filetype
  data, ndata = readbinary(filename, isprint=isprint)

  if ndata != int(np.prod(shape)):
    err = str(ndata)+" is incorrect data length for "+var+\
          " defined for "+str(ntime)+" timesteps"+\
         " on grid with dims = "+str(ndims)
    raise ValueError(err)
  else:
    data = np.reshape(np.asarray(data, dtype=dtype), shape)
  
  return data

def read_dimless_thermodynamics_binary(thermofile, ndims,
                                       ntime, nspacedims):

  # expected lengths of data defined on gridbox centres or faces 
  cen = [ntime, int(np.prod(ndims))]
  zface = [ntime, int(ndims[2] * ndims[1] * (ndims[0]+1))]
  xface = [ntime, int(ndims[2] * (ndims[1]+1) * ndims[0])]
  yface = [ntime, int((ndims[2]+1) * ndims[1] * ndims[0])]

  thermodata = {} 
  
  vars = ["press", "temp", "qvap", "qcond"]
  datatypes = [np.double]*4
  for v, var in enumerate(vars): 
    thermodata[var] = thermovar_from_binary(var, thermofile, cen,
                                            ntime, ndims, datatypes[v],
                                            isprint=False)

  datatypes = [np.double]*3
  if nspacedims >= 1:
    thermodata["wvel"] = thermovar_from_binary("wvel", thermofile, zface,
                                              ntime, ndims, datatypes[0],
                                              isprint=False)
    if nspacedims >= 2:
      thermodata["uvel"] = thermovar_from_binary("uvel", thermofile, xface,
                                                ntime, ndims, datatypes[1],
                                                isprint=False)
      if nspacedims >= 3:
        thermodata["vvel"] = thermovar_from_binary("vvel", thermofile, yface,
                                                  ntime, ndims, datatypes[2],
                                                  isprint=False)

  return thermodata

def get_thermodynamics_from_thermofile(thermofile, ndims, inputs=False,
                                       constsfile="", configfile=""):

  if not inputs:  
    inputs = thermoinputsdict(configfile, constsfile)

  thermodata = read_dimless_thermodynamics_binary(thermofile, ndims, 
                                                  inputs["ntime"],
                                                  inputs["nspacedims"]) # dimensionless data [time, gridboxes]
  thermodata = ThermoOnGrid(thermodata, inputs, ndims) 

  return thermodata # data with units in 4D arrays with dims [time, y, x, z]

def plot_thermodynamics(constsfile, configfile, gridfile,
                                  thermofile, binpath, savefig):

    plt.rcParams.update({'font.size': 14})

    inputs = thermoinputsdict(configfile, constsfile)
    gbxbounds, ndims = rgrid.read_dimless_gbxboundaries_binary(gridfile,
                                                            COORD0=inputs["COORD0"],
                                                            return_ndims=True,
                                                            isprint=False)
    xyzhalf = rgrid.halfcoords_from_gbxbounds(gbxbounds, isprint=False) #[m]
    zhalf, xhalf, yhalf = [half/1000 for half in xyzhalf] #convery [m] to [km]
    zfull, xfull, yfull = rgrid.fullcell_fromhalfcoords(zhalf, xhalf, yhalf) #[m]

    thermodata = get_thermodynamics_from_thermofile(thermofile, ndims,
                                                    inputs=inputs)

    plot_1dprofiles(zfull, thermodata, inputs["Mr_ratio"], binpath, savefig)

    if inputs["nspacedims"] > 1: 
      xxh, zzh = np.meshgrid(xhalf, zhalf, indexing="ij") # dims [xdims, zdims]
      xxf, zzf = np.meshgrid(xfull, zfull, indexing="ij") # dims [xdims, zdims]
      plot_2dcolormaps(zzh, xxh, zzf, xxf, thermodata, inputs, binpath, savefig)
      plot_2dwindfield(zzh, xxh, zzf, xxf, thermodata["wvel_cens"],
                      thermodata["uvel_cens"], binpath, savefig)

def try1dplot(ax, data, zfull):
  try:
    ax.plot(data, zfull, marker="x") # (fails for 0D model)
  except:
    ax.scatter(data, zfull, marker="x")

def plot_1dprofiles(zfull, thermodata, Mr_ratio, binpath, savefig):
    
    vars = ["press", "temp", "qvap", "qcond",
            "wvel_cens", "uvel_cens", "vvel_cens"]
    units = [" /Pa", " /K", "", ""]+[" /ms$^{-1}$"]*3
    
    fig, axs = plt.subplots(nrows=2, ncols=4, figsize=(16, 8))
    axs = axs.flatten()

    lines = 0
    for v, var in enumerate(vars):
      if var in thermodata.vars:
        profalltime = thermodata.xymean(thermodata[var]) # 1d profile at all times
        try1dplot(axs[v], profalltime.T, zfull[None,:].T)
        axs[v].set_xlabel(vars[v]+units[v])  
        lines += 1

    supersat = relative_humidity(thermodata.xymean(thermodata.press),
                                  thermodata.xymean(thermodata.temp),
                                  thermodata.xymean(thermodata.qvap),
                                  Mr_ratio)[1]
    try1dplot(axs[lines], supersat.T, zfull[None,:].T)
    axs[lines].set_xlabel("supersaturation")  
    lines+=1

    for a in range(lines, len(axs), 1):
      axs[a].remove() # delete unused axes
    axs[0].set_ylabel("z /km")
    axs[4].set_ylabel("z /km")

    fig.tight_layout()
    if savefig:
        savename = "/thermo1dalltimeprofiles.png"
        fig.savefig(binpath+savename, dpi=400,
                    bbox_inches="tight", facecolor='w', format="png")
        print("Figure .png saved as: "+binpath+savename)
    plt.show()

def plot_2dcolormaps(zzh, xxh, zzf, xxf,
                     thermodata, inputs, binpath, savefig):

  vars = ["press", "temp", "qvap", "qcond"]
  units = [" /Pa", " /K", "", ""]
  cmaps = ["PRGn", "RdBu_r", "BrBG", "BrBG"]
  cmapcens = [np.nanmin(thermodata.press),
              np.nanmin(thermodata.temp),
              np.mean(thermodata.qvap), 0.0] 

  fig, axs = plt.subplots(nrows=2, ncols=3,  figsize=(13,8))
  axs = axs.flatten()
  
  for v, var in enumerate(vars):
    mean2d = thermodata.ytmean(thermodata[var]) #avg over time and y axes
    norm=colors.CenteredNorm(vcenter=cmapcens[v])
    pcm = axs[v].pcolormesh(xxh[:,:], zzh[:,:], mean2d, cmap=cmaps[v], norm=norm)
    plt.colorbar(pcm, ax=axs[v], location="top", label=var+units[v])

  relh, supersat = relative_humidity(thermodata.press, thermodata.temp,
                                      thermodata.qvap, inputs["Mr_ratio"])
  
  cmaps = ["Blues", "PuOr"]
  relh_supersat_colomaps([axs[4], axs[5]], zzh, xxh, zzf, xxf, cmaps,
                          thermodata.ytmean(relh), thermodata.ytmean(supersat))
  
  for ax in axs:
    ax.set_aspect("equal")

  axs[0].set_ylabel("z /km")
  axs[3].set_ylabel("z /km")
  for ax in axs[3:]:
    ax.set_xlabel("x /km")

  fig.tight_layout()
  if savefig:
      savename = "/thermo2dmeanprofiles.png"
      fig.savefig(binpath+savename, dpi=400,
                  bbox_inches="tight", facecolor='w', format="png")
      print("Figure .png saved as: "+binpath+savename)
  plt.show()

def relh_supersat_colomaps(axs, zzh, xxh, zzf, xxf, cmaps,
                           relh, supersat):
  
  relh = relh * 100 # convert relative humidity to %
  
  label = ["% relative humidity", "supersaturation"]
  norms = [colors.CenteredNorm(vcenter=np.mean(relh)),
            colors.TwoSlopeNorm(vcenter=0.0)] 
  contour = [1.0, 0.0]
  for v, var in enumerate([relh, supersat]):
    pcm = axs[v].pcolormesh(xxh, zzh, var, cmap=cmaps[v], norm=norms[v])                  
    cb = plt.colorbar(pcm, ax=axs[v], location="top", label=label[v])

    cbticks = np.linspace(pcm.norm.vmin, pcm.norm.vmax, 5)
    if v==1: cbticks = [pcm.norm.vmin, pcm.norm.vcenter, pcm.norm.vmax]
    cb.ax.set_xticks(cbticks) 
    cb.ax.set_xticklabels(["{:.3g}".format(c) for c in cbticks]) 
    
    axs[v].contour(xxf, zzf, var, levels=[contour[v]],
                  linestyles=["--"], colors=["grey"])
    cb.ax.plot([contour]*2, [0, 1], color='grey', linewidth=0.95)

def plot_2dwindfield(zzh, xxh, zzf, xxf, wvel_cens, uvel_cens, 
                     binpath, savefig):
  
  wcen = np.mean(wvel_cens, axis=(0,1)) # avg over y and time axes
  ucen = np.mean(uvel_cens, axis=(0,1))
  norm = np.sqrt(wcen**2 + ucen**2)

  fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(9,6))
  axs = axs.flatten()

  label = ["w", "u", "|wind velocity|"]
  units = " / m s$^{-1}$"
  cmaps = ["coolwarm"]*3
  for v, vel in enumerate([wcen, ucen, norm]):
    pcm = axs[v].pcolormesh(xxh, zzh, vel, cmap=cmaps[v])
    plt.colorbar(pcm, ax=axs[v], location="top", label=label[v]+units) 
  axs[2].quiver(xxf, zzf, ucen, wcen)
  
  axs[0].set_ylabel("z /km")
  for ax in axs:
    ax.set_xlabel("x /km")
    ax.set_aspect("equal")

  fig.tight_layout()
  if savefig:
      savename = "/thermowindprofiles.png"
      fig.savefig(binpath+savename, dpi=400,
                  bbox_inches="tight", facecolor='w', format="png")
      print("Figure .png saved as: "+binpath+savename)
  plt.show()