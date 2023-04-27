import numpy as np
import matplotlib.pyplot as plt

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

class Thermo4D:

    def __init__(self, thermodata, ndims, ntime):
        
        reshape = [ntime] + list(np.flip(ndims))
        
        self.vars = ["press", "temp", "qvap", "qcond"]
        self.press = np.reshape(thermodata["press"], reshape)
        self.temp = np.reshape(thermodata["temp"], reshape)
        self.qvap= np.reshape(thermodata["qvap"], reshape)
        self.qcond = np.reshape(thermodata["qcond"], reshape)
        # self.wvel = np.array([])
        # self.uvel = np.array([])
        # self.vvel = np.array([])

        if "wvel" in thermodata.keys():
            self.wvel = np.reshape(thermodata["wvel"], reshape) 
            self.vars += ["wvel"]
            if "uvel" in thermodata.keys():
                self.uvel = np.reshape(thermodata["uvel"], reshape) 
                self.vars += ["uvel"]
                if "vvel" in thermodata.keys():
                    self.vvel = np.reshape(thermodata["vvel"], reshape)  
                    self.vars += ["vvel"]
    
    def __getitem__(self, key):
      if key not in self.vars:
        err = "no variable "+key+" in Thermo4D"
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

    def meanytime(self, var):

        return np.mean(var, axis=(0,1))

    def meanxytime(self, var):

        return np.mean(var, axis=(0,1,2))

def thermovar_from_binary(var, thermofiles, shape,
                          ntime, ndims, dtype):
  
  filestem, filetype = thermofiles.split(".")
  filename = filestem+"_"+var+"."+filetype
  data, ndata = readbinary(filename)

  if ndata != int(np.prod(shape)):
    err = "incorrect data length for "+var+" defined at"+\
                    " gridbox centres for "+str(ntime)+" timesteps"+\
                    " on grid with dims = "+str(ndims)
    raise ValueError(err)
  else:
    data = np.reshape(np.asarray(data, dtype=dtype), shape)
  
  return data

def read_dimless_thermodynamics_binary(thermofiles, ndims,
                                       ntime, SDnspace):

  # expected lengths of data defined on gridbox centres or faces 
  cen = [ntime, int(np.prod(ndims))]
  zface = [ntime, int((ndims[0]+1)*ndims[1]*ndims[2])]
  xface = [ntime, int((ndims[1]+1)*ndims[2]*ndims[0])]
  yface = [ntime, int((ndims[2]+1)*ndims[0]*ndims[1])]

  thermodata = {} 
  
  vars = ["press", "temp", "qvap", "qcond"]
  datatypes = [np.double]*4
  for v, var in enumerate(vars): 
    thermodata[var] = thermovar_from_binary(var, thermofiles, cen,
                                            ntime, ndims, datatypes[v])

  datatypes = [np.double]*3
  if SDnspace >= 1:
    thermodata["wvel"] = thermovar_from_binary("wvel", thermofiles, zface,
                                              ntime, ndims, datatypes[0]) 
    if SDnspace >= 2:
      thermodata["uvel"] = thermovar_from_binary("uvel", thermofiles, xface,
                                                ntime, ndims, datatypes[1])
      if SDnspace >= 3:
        thermodata["vvel"] = thermovar_from_binary("vvel", thermofiles, yface,
                                                  ntime, ndims, datatypes[2])

  return thermodata

def get_thermodynamics_from_thermofile(thermofile, ndims, inputs=False,
                                       constsfile="", configfile=""):

  if not inputs:  
    inputs = thermoinputsdict(configfile, constsfile)

  thermodata = read_dimless_thermodynamics_binary(thermofile, ndims, 
                                                  inputs["ntime"],
                                                  inputs["SDnspace"]) 
  dth = DimlessThermodynamics(inputs=inputs)
  thermodata = dth.redimensionalise(thermodata)

  thermodata = Thermo4D(thermodata, ndims, inputs["ntime"]) 

  return thermodata # shaped into 4D array dims [time, y, x, z]

def plot_thermodynamics_timeslice(constsfile, configfile, gridfile,
                                  thermofile, binpath, savefig):

    plt.rcParams.update({'font.size': 14})

    inputs = thermoinputsdict(configfile, constsfile)
    gbxbounds, ndims = rgrid.read_dimless_gbxboundaries_binary(gridfile,
                                                            COORD0=inputs["COORD0"],
                                                            return_ndims=True)
    zhalf, xhalf, yhalf = rgrid.halfcoords_from_gbxbounds(gbxbounds)
    zfull, xfull, yfull = rgrid.fullcell_fromhalfcoords(zhalf, xhalf, yhalf)

    thermodata = get_thermodynamics_from_thermofile(thermofile, ndims,
                                                    inputs=inputs)
    
    plot_1dprofiles(zfull, thermodata, binpath, savefig)
    
    xxh, zzh = np.meshgrid(xhalf, zhalf, indexing="ij") # dims [xdims, zdims]
    xxf, zzf = np.meshgrid(xfull, zfull, indexing="ij") # dims [xdims, zdims]
    plot_2dcolormaps(zzh, xxh, zzf, xxf, thermodata, inputs, binpath, savefig)
    plot_2dwindfield(zzh, xxh, zzf, xxf, thermodata["wvel"],
                    thermodata["uvel"], binpath, savefig)

def plot_1dprofiles(zfull, thermodata, binpath, savefig):
    
    vars = thermodata.vars
    units = [" /Pa", " /K", "", ""]+[" /ms$^{-1}$"]*3
    fig, axs = plt.subplots(nrows=2, ncols=4, figsize=(12, 6))
    axs = axs.flatten()

    handles, labels = [], []
    for v, var in enumerate(vars):
      var1d_alltime = np.mean(thermodata[var], axis=(1,2)) # z profile at each time
      line = axs[v].plot(var1d_alltime.T, zfull[None,:].T, marker="x")
      axs[v].set_xlabel(vars[v]+units[v])  
      handles.append(line)

    axs[0].set_ylabel("z /km")
    axs[4].set_ylabel("z /km")
    for a in range(len(vars), len(axs), 1):
      axs[a].remove()

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
  cmaps = ["PRGn", "RdBu", "BrBG", "BrBG", "Blues", "Blues"]
  
  fig, axs = plt.subplots(nrows=2, ncols=3,  figsize=(8,6))
  axs = axs.flatten()
  
  for v, var in enumerate(vars):
    mean2d = thermodata.meanytime(thermodata[var]) #avg over time and y axes
    pcm = axs[v].pcolormesh(xxh[:,:], zzh[:,:], mean2d, cmap=cmaps[v])
    plt.colorbar(pcm, ax=axs[v], location="top", label=var+units[v])

  relh, supersat = relative_humidity(thermodata.press, thermodata.temp,
                                      thermodata.qvap, inputs["Mr_ratio"])
  
  meanrelh = thermodata.meanytime(relh)
  pcm =axs[4].pcolormesh(xxh[:,:], zzh[:,:], meanrelh, cmap=cmaps[4])                  
  cb = plt.colorbar(pcm, ax=axs[4], location="top", label="relative humidity")
  axs[4].contour(xxf, zzf, meanrelh, levels=[1.0],
                linestyles=["--"], colors=["grey"])
  cb.ax.plot([1.0]*2, [0, 1], color='grey', linewidth=0.95)

  meansupersat = thermodata.meanytime(supersat)*100
  pcm = axs[5].pcolormesh(xxh[:,:], zzh[:,:], meansupersat, cmap=cmaps[5])
  cb = plt.colorbar(pcm, ax=axs[5], location="top", label="% supersaturation")
  axs[5].contour(xxf, zzf, meansupersat, levels=[0.0],
                linestyles=["--"], colors=["grey"])
  cb.ax.plot([0.0]*2, [0, 1], color='grey', linewidth=0.95)

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

def plot_2dwindfield(zzh, xxh, zzf, xxf, wvel4d, uvel4d, 
                     binpath, savefig):
  
  cmaps = ["coolwarm"]*3
  meanwvel = np.mean(wvel4d, axis=(0,1)) #avg over y and time axes
  meanuvel = np.mean(uvel4d, axis=(0,1))
  norm = np.sqrt(meanwvel**2 + meanuvel**2)

  fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(9,6))
  axs = axs.flatten()

  pcm = axs[0].pcolormesh(xxh[:,:], zzh[:,:], meanwvel, cmap=cmaps[0])
  plt.colorbar(pcm, ax=axs[0], location="top", label="w / m s$^{-1}$") 

  pcm = axs[1].pcolormesh(xxh[:,:], zzh[:,:], meanuvel, cmap=cmaps[0])
  plt.colorbar(pcm, ax=axs[1], location="top", label="u / m s$^{-1}$")

  pcm = axs[2].pcolormesh(xxh[:,:], zzh[:,:], norm, cmap=cmaps[0])
  plt.colorbar(pcm, ax=axs[2], location="top",
               label="|wind velocity| / m s$^{-1}$")
  axs[2].quiver(xxf, zzf, meanuvel, meanwvel)
  
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