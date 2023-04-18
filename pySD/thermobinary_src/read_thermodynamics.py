import numpy as np
import matplotlib.pyplot as plt

from .create_thermodynamics import thermoinputsdict, DimlessThermodynamics
from ..gbxboundariesbinary_src.read_gbxboundaries import get_fullcoords_from_gridfile
from ..readbinary import readbinary

def read_dimless_thermodynamics_binary(thermofile, ngridboxes, ntime):

  datatypes = [np.double]*7

  thermodata = {
    "press": [],
    "temp": [],
    "qvap": [],
    "qcond": [],
    "wvel": [],
    "vvel": [],
    "uvel": []
  }
  
  filestem, filetype = thermofile.split(".")
  for v, var in enumerate(thermodata.keys()):
    filename = filestem+"_"+var+"."+filetype
    data, ndata = readbinary(filename)

    if ndata != [ngridboxes*ntime]:
      raise ValueError("incorrect data length for "+var+" given "+
                       str(ngridboxes)+ " gridboxes and "+str(ntime)+\
                        " timesteps")
    data = np.asarray(data, dtype=datatypes[v])
    thermodata[var] = np.reshape(data, [ngridboxes, ntime])
  
  return thermodata


def get_thermodyanmcis_from_thermofile(thermofile, ngridboxes,
                                       constsfile, configfile):
    
  inputs = thermoinputsdict(configfile, constsfile)

  
  thermodata = read_dimless_thermodynamics_binary(thermofile, ngridboxes, 
                                                  inputs["ntime"]) 
  dth = DimlessThermodynamics(inputs=inputs)
  
  return dth.redimensionalise(thermodata)

def plot_thermodynamics_timeslice(constsfile, configfile, gridfile,
                                  thermofile, t2plt, binpath, savefig):

    plt.rcParams.update({'font.size': 14})

    zfull, xfull, yfull = get_fullcoords_from_gridfile(gridfile,
                                                       constsfile=constsfile)
    
    ngridboxes = len(zfull)*len(xfull)*len(yfull)
    thermodata = get_thermodyanmcis_from_thermofile(thermofile, ngridboxes,
                                                    constsfile, configfile)
    
    vars = ["press", "temp", "qvap", "qcond", "wvel", "vvel", "uvel"]
    units = [" /Pa", " /K", "", ""]+[" /ms$^{-1}$"]*3
    fig, axs = plt.subplots(nrows=2, ncols=4, figsize=(12, 6))
    axs = axs.flatten()

    for v, var in enumerate(vars):
      axs[v].plot(thermodata[var], range(ngridboxes))
      axs[v].set_xlabel(var+units[v])
      axs[v].set_ylabel("GBx")
    axs[v+1].remove()

    fig.tight_layout()
    if savefig:
        tstr = str(t2plt).replace('.', "p")
        fig.savefig(binpath+"/thermodynamics_time"+tstr+".png", dpi=400,
                    bbox_inches="tight", facecolor='w', format="png")
        print("Figure .png saved as: "+binpath+"/gridboxboundaries.png")
    plt.show()
