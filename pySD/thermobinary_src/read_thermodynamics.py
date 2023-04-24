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

def read_dimless_thermodynamics_binary(thermofile, ngridboxes,
                                       ntime, SDnspace):

  datatypes = [np.double]*7
  vars = ["press", "temp", "qvap", "qcond"]
  if SDnspace >= 1:
    vars += ["wvel"]
    if SDnspace >= 2:
      vars += ["uvel"]
      if SDnspace >= 3:
        vars += ["vvel"]

  thermodata = {} 
  filestem, filetype = thermofile.split(".")
  for v, var in enumerate(vars):
    filename = filestem+"_"+var+"."+filetype
    data, ndata = readbinary(filename)

    if ndata != [ngridboxes*ntime]:
      raise ValueError("incorrect data length for "+var+" given "+
                      str(ngridboxes)+ " gridboxes and "+str(ntime)+\
                        " timesteps")
    data = np.asarray(data, dtype=datatypes[v])
    thermodata[var] = np.reshape(data, [ntime, ngridboxes])
  
  return thermodata

def get_thermodynamics_from_thermofile(thermofile, ngridboxes, inputs=False,
                                       constsfile="", configfile=""):

  if not inputs:  
    inputs = thermoinputsdict(configfile, constsfile)

  thermodata = read_dimless_thermodynamics_binary(thermofile, ngridboxes, 
                                                  inputs["ntime"],
                                                  inputs["SDnspace"]) 
  dth = DimlessThermodynamics(inputs=inputs)
  
  return dth.redimensionalise(thermodata)

def plot_thermodynamics_timeslice(constsfile, configfile, gridfile,
                                  thermofile, t2plt, binpath, savefig):

    plt.rcParams.update({'font.size': 14})

    inputs = thermoinputsdict(configfile, constsfile)
    gbxbounds, ndims = rgrid.read_dimless_gbxboundaries_binary(gridfile,
                                                            COORD0=inputs["COORD0"],
                                                            return_ndims=True)
    zhalf, xhalf, yhalf = rgrid.halfcoords_from_gbxbounds(gbxbounds)
    zfull, xfull, yfull = rgrid.fullcell_fromhalfcoords(zhalf, xhalf, yhalf)

    ngridboxes = int(np.prod(ndims))
    thermodata = get_thermodynamics_from_thermofile(thermofile, ngridboxes,
                                                    inputs=inputs)
    
    vars = [v for v in thermodata.keys()]
    units = [" /Pa", " /K", "", ""]+[" /ms$^{-1}$"]*3
    fig, axs = plt.subplots(nrows=2, ncols=4, figsize=(12, 6))
    axs = axs.flatten()

    gbxs = range(ngridboxes)
    if t2plt == "all":
      for v, var in enumerate(vars):
        for t in range(inputs['ntime']):
          axs[v].plot(thermodata[var][t,:], gbxs, marker="x")
      
    else:
      for v, var in enumerate(vars):
        end, stp = inputs["T_END"]+inputs["COUPLTSTEP"], inputs["COUPLTSTEP"]
        times = np.arange(0, end, stp)
        t = np.argmin(abs(t2plt-times))
        l = axs[v].plot(thermodata[var][t,:], gbxs,
                        color="k", marker="x", linestyle="--")
      
        axs[0].legend(["t={:.0f}s".format(times[t])], loc="upper right")
    
    for v in range(len(vars)):
      axs[v].set_xlabel(vars[v]+units[v])
      axs[v].set_ylabel("GBx")
    for a in range(len(vars), len(axs), 1):
      axs[a].remove()

    fig.tight_layout()
    if savefig:
        tstr = str(t2plt).replace('.', "p")
        fig.savefig(binpath+"/thermodynamics_time"+tstr+".png", dpi=400,
                    bbox_inches="tight", facecolor='w', format="png")
        print("Figure .png saved as: "+binpath+"/gridboxboundaries.png")
    plt.show()
