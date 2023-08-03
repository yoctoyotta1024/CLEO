### Author: Clara Bayley
### File: exmpl_figsrc.py
### classes and functions called by exmpl_plot.py
### to plot figures of some of the output
### data from the example run of CLEO

import sys
import numpy as np
import xarray as xr
import awkward as ak
import random
import matplotlib.pyplot as plt
from matplotlib.markers import MarkerStyle

sys.path.append("..") # used to include directory containing pySD package in python interpreter PATH
from pySD import cxx2py
import pySD.gbxboundariesbinary_src.read_gbxboundaries as rgrid

### ------------------- post-processing functions ------------------ ###

def get_rawdataset(dataset):
  print("Dataset: "+dataset)
  return xr.open_dataset(dataset, engine="zarr", consolidated=False)

def raggedvar_fromzarr(ds, raggedcount, var):
  ''' returns ragged ak.Array dims [time, ragged]
  for a variable "var" in zarr ds '''
  return ak.unflatten(ds[var].values, raggedcount)

def var4d_fromzarr(ds, ntime, ndims, key):
  '''' returns 4D variable with dims
  [time, y, x, z] from zarr dataset "ds" '''
  
  reshape = [ntime] + list(ndims)
  return np.reshape(ds[key].values, reshape) 

def get_timemins(dataset):
  ''' returns time in minutes from
  dataset containing time in seconds'''

  time = get_rawdataset(dataset)["time"].values # assumed to be [seconds]
  return time / 60 # convert seconds to minutes

class Setup:
  ''' class for handling dictionary
  containing configuration from setuptxt '''

  def __init__(self, setuptxt, ntime=None, isprint=True):

    print("Reading setuptxt file:\n  "+setuptxt)
    self.setupdict = cxx2py.read_configtxt_into_floats(setuptxt, False)[0]
    self.setupdict.update(cxx2py.read_cpp_into_floats(setuptxt, False)[0])
    self.setupdict.update(cxx2py.derive_more_floats(self.setupdict, False))
   
    if not ntime:
      ntime = int(np.ceil(self.setupdict["T_END"]/self.setupdict["OBSTSTEP"]))
    self.setupdict["ntime"] = ntime
    
    if isprint:
      cxx2py.print_dict_statement(setuptxt, "setup", self.setupdict)

  def __getitem__(self, key):
    return self.setupdict[key]

class GridBoxes:
  ''' domain info and (x, y, z) coords for gridboxes
  as well as volume info and 2D (z,x) meshgrids and areas '''

  def __init__(self, gridfile, setup, isprint=True):
    
    COORD0 = setup["COORD0"]
    gbxbounds, ndims = rgrid.read_dimless_gbxboundaries_binary(gridfile,
                                                              COORD0=COORD0,
                                                              return_ndims=True,
                                                              isprint=isprint) 
    domainvol, gbxvols, ngrid = rgrid.domaininfo(gbxbounds,
                                                 isprint=isprint)
 
    self.ngrid = ngrid          # number of gridboxes in domain
    self.ndims = np.flip(ndims) # domain dimensions (ie. no. gridboxes in [y,x,z] direction)
    
    halfcoords = rgrid.halfcoords_from_gbxbounds(gbxbounds, isprint=False)
    self.zhalf, self.xhalf, self.yhalf = halfcoords
    self.zfull = rgrid.fullcell(self.zhalf)  # full cell z coords (centres)
    self.xfull = rgrid.fullcell(self.xhalf)  # full cell x coords (centres)
    self.yfull =  rgrid.fullcell(self.yhalf) # full cell y coords (centres)

    self.domainvol = domainvol                     # total volume of domain
    self.gbxvols = np.reshape(gbxvols, self.ndims) # volume of each gridbox in domain
    self.gbxareas = self.gbxvols / np.diff(self.zhalf)[None, None, :] # x-y plane horizontal areas

    self.xxh, self.zzh = np.meshgrid(self.xhalf, self.zhalf, indexing="ij") # dims [xhalf, zhalf] [m]
    self.xxf, self.zzf = np.meshgrid(self.xfull, self.zfull, indexing="ij") # dims [xfull, zfull] [m]

class MassMoments:

  def __init__(self, dataset, setup, ndims):
    
    ds =  get_rawdataset(dataset)
    ntime = setup["ntime"]
    
    self.nsupers =  var4d_fromzarr(ds, ntime, ndims, "nsupers")        # number of superdroplets in gbxs over time
    self.mom0 = var4d_fromzarr(ds, ntime, ndims, "mom0")               # number of droplets in gbxs over time
    self.mom1 = var4d_fromzarr(ds, ntime, ndims, "mom1")               # total mass of droplets in gbxs over time
    self.mom2 = var4d_fromzarr(ds, ntime, ndims, "mom2")               # 2nd mass moment of droplets (~reflectivity)
    
    divmom1 = np.where(self.mom1 > 0, self.mom1, np.inf)
    self.effmass = self.mom2 / divmom1                               # Effective Radius of droplets

class SDData:
  
  def __init__(self, dataset):
    
    ds = get_rawdataset(dataset)
     
    self.totnsupers = ds["raggedcount"].values # total no. SDs in domain at each tstep

    self.sdindex = raggedvar_fromzarr(ds, self.totnsupers, "sdindex")
    self.eps = raggedvar_fromzarr(ds, self.totnsupers, "eps")
    self.radius = raggedvar_fromzarr(ds, self.totnsupers, "radius")
    self.m_sol = raggedvar_fromzarr(ds, self.totnsupers, "m_sol")

    self.coord3 = raggedvar_fromzarr(ds, self.totnsupers, "coord3")
    self.coord1 = raggedvar_fromzarr(ds, self.totnsupers, "coord1")
  
  def __getitem__(self, key):
    if key == "totnsupers":
      return self.totnsupers
    elif key == "sdindex":
      return self.sdindex
    elif key == "eps":
      return self.eps
    elif key == "radius":
      return self.radius
    elif key == "m_sol":
      return self.m_sol
    elif key == "coord3":
      return self.coord3
    elif key == "coord1":
      return self.coord1

class SurfPrecip:

  def __init__(self, dataset, setup, gbxs):
    
    ds = get_rawdataset(dataset)
    surfpp = ds["surfpp"].values # [mm]
    deltat = np.mean(np.diff(ds["time"].values)) / 60 / 60 # [hrs]
    
    reshape = [setup["ntime"]] + list(gbxs.ndims)[0:2]
    self.surfpp = np.reshape(surfpp, reshape) # dims [time, x, y]
    
    self.accum = np.cumsum(self.surfpp, axis=0) # [mm]
    self.rate = self.surfpp / deltat # [mm/hr]

    gbxareas = gbxs.gbxareas[:,:,0] # areas of surface gbxs
    scale_area_factor = gbxareas / np.sum(gbxareas) # total area of surface ( = gbxs["domainarea"])
    self.totrate = np.sum(self.rate * scale_area_factor, axis=(1,2))
    self.totaccum = np.sum(self.accum * scale_area_factor, axis=(1,2)) 

### ---------------------------------------------------------------- ###

### ---------------------- plotting functions ---------------------- ###

def save_figure(fig, savedir, savename, show=True):

  fig.savefig(savedir+savename, 
              dpi=400, bbox_inches="tight", 
              facecolor='w', format="png")
  print("Figure .png saved as: "+savedir+savename)

  if show:
      plt.show()

def genericplot(plot_function, argsdict, figsize,
                savefig=False, savedir=None, savename=None):
  ''' shows and optionally saves a figure of figsize
  by calling the plot_function with the arguments given
  by argsdict. Given that signature of plot_funtion is
  plot_function(figsize, {args}), return [fig, axs, lines]'''

  fig, axs, lines = plot_function(figsize, **argsdict)

  fig.tight_layout()
  if savefig:
    save_figure(fig, savedir, savename, show=True)
  else:
    plt.show()

  return fig, axs, lines

def plot_domainmassmoments(figsize, time, massmoms):
 
  def totmassmom(massmom):
    '''mass moment summed over entire domain'''
    return  np.sum(massmom, axis=(1,2,3))
  
  fig, axs = plt.subplots(nrows=2, ncols=2, figsize=figsize,
                          sharex=True)
  axs = axs.flatten() 
  fig.suptitle("Mass Moments Over Entire Domain")
  
  l0 = axs[0].plot(time, totmassmom(massmoms.mom0))
  l1 = axs[1].plot(time, totmassmom(massmoms.mom1)) 
  l2 = axs[2].plot(time, totmassmom(massmoms.mom2))
  meaneffmass = np.mean((massmoms.effmass), axis=(1,2,3))
  l3 = axs[3].plot(time, meaneffmass)

  axs[0].set_ylabel("$\u03BB^{m}_{0}$, number of  droplets")
  axs[1].set_ylabel("$\u03BB^{m}_{1}$, droplet mass /g")
  axs[2].set_ylabel("$\u03BB^{m}_{2}$ ~reflectivity /g$^2$")
  ylab3 = "mean effective droplet mass,\n<$\u03BB^{m}_{2}$/$\u03BB^{m}_{1}>$ /g"
  axs[3].set_ylabel(ylab3)
  
  for ax in axs:
    ax.set_xlabel("time /min")
  
  return fig, axs, [l0, l1, l2, l3]

def attrtimeseries_for_ids(sddata, attr, ids):
  ''' returns 2D array with dimensions [time, SD]
  containing attribute data over time for 'ndrops'
  randomly selected from superdrops with id in
  range [minid, maxid] '''
    
  def attrtimeseries(sddata, attr, id):
    '''selects attribute from sddata belonging
    to superdroplet with identitiy 'sdindex'
    at every output time '''

    bools = ak.Array(sddata.sdindex==id) # True/False id is found in sdindex at each time
    attr_timeseries = sddata[attr][bools]
    num = ak.num(attr_timeseries) # at each time, number of positions where id is found (should be 0 or 1)
    
    if any(num[num!=1]):
      errmsg = "attr_timeseries has times when more"+\
        " than one position in sddata has sdindex==id. num should be"+\
      " list of either 1 or 0 (id found in sdindex at given time or not)"
      raise ValueError(errmsg)
    else:
      attr_timeseries = np.where(num!=0, attr_timeseries, ak.Array([[np.nan]])) # replace empty positions with nan
      return ak.flatten(attr_timeseries, axis=1)
  
  ndrops_attr = []
  for id in ids: 
    ndrops_attr.append(attrtimeseries(sddata, attr, id))
  
  return np.asarray(ndrops_attr).T

def plot_randomsampleSDs(figsize, time, sddata, nsample):

  fig, axs = plt.subplots(nrows=2, ncols=3, figsize=figsize)    
  axs = axs.flatten() 
  fig.suptitle("Time Evolution of a Random Sample of Superdroplets")

  minid, maxid = 0, sddata.totnsupers[0] # largest value of ids to sample
  ids2plot = random.sample(list(range(minid, maxid, 1)), nsample)
  
  radii = attrtimeseries_for_ids(sddata, "radius", ids=ids2plot) 
  m_sols = attrtimeseries_for_ids(sddata, "m_sol", ids=ids2plot)
  epss = attrtimeseries_for_ids(sddata, "eps", ids=ids2plot)
  zkms = attrtimeseries_for_ids(sddata, "coord3", ids=ids2plot) / 1000 # [km] 
  xkms = attrtimeseries_for_ids(sddata, "coord1", ids=ids2plot) / 1000 # [km] 

  l0s = axs[0].plot(time, radii, linewidth=0.8)
  l1s = axs[1].plot(time, epss, linewidth=0.8)
  l2s = axs[2].plot(time, m_sols, linewidth=0.8)
  mks = MarkerStyle('o', fillstyle='full')
  l3s = axs[3].plot(time, zkms, linestyle="", marker=mks, markersize=0.2)
  l4s = axs[4].plot(time, xkms, linestyle="", marker=mks, markersize=0.2)
  l5s = axs[5].plot(xkms, zkms, linestyle="", marker=mks, markersize=0.2) 

  axs[0].set_yscale('log')
  axs[0].set_ylabel('radius /\u03BCm')
  axs[1].set_ylabel('multiplicity, \u03BE')
  axs[2].set_ylabel('solute mass /g')
  axs[3].set_ylabel('zcoord /km')
  axs[4].set_ylabel('xcoord /km')
  for ax in axs[0:5]:
    ax.set_xlabel("time /min")
  
  axs[5].set_ylabel('zcoord /km')
  axs[5].set_xlabel("xcoord /km")
  axs[5].set_aspect("equal")

  return fig, axs, [l0s, l1s, l2s, l3s, l4s, l5s] 

def plot_surfaceprecip(figsize, time, precip):

  fig, axs = plt.subplots(nrows=1, ncols=2,
                          figsize=figsize, sharex=True)
  axs = axs.flatten()
  fig.suptitle("Surface Precipitation over Entire Domain")

  l0 = axs[0].plot(time, precip.totrate, color="k")
  l1 = axs[1].plot(time, precip.totaccum, color="k")

  axs[0].set_ylabel("rate / mm hr$^{-1}$")
  axs[1].set_ylabel("accumulated /mm")
  for ax in axs:
    ax.set_xlabel("time /min")

  return fig, axs, [l0, l1]
### ---------------------------------------------------------------- ###