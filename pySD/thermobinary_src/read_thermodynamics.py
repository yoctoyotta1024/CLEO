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
    
    print(thermodata)
    fig, axs = plt.subplots(nrows=3, ncols=4, figsize=(10, 5))

    # for i, crd in enumerate(["z", "x", "y"]):
    #     xlims = [0.8*np.amin(deltas[i])/1000, np.amax(deltas[i])/1000*1.2]
    #     ylims = [np.amin(halfs[i])/1000, np.amax(halfs[i])/1000]
    #     axs[i].scatter(deltas[i]/1000, fulls[i]/1000,
    #                    color="k", label="centres")
    #     axs[i].hlines(halfs[i]/1000, xlims[0], xlims[1],
    #                   color="grey", alpha=0.8, linewidth=0.8, label="boundaries")
    #     axs[i].set_xlim(xlims)
    #     axs[i].set_ylim(ylims)
    #     axs[i].set_xlabel("gridbox spacing, \u0394 "+crd+" /km")
    #     axs[i].set_ylabel("gridbox centres, "+crd+"f /km")
    #     axs[i].legend()

    # fig.tight_layout()
    # if savefig:
    #     fig.savefig(binpath+"/gridboxboundaries.png", dpi=400,
    #                 bbox_inches="tight", facecolor='w', format="png")
    #     print("Figure .png saved as: "+binpath+"/gridboxboundaries.png")
    # plt.show()
