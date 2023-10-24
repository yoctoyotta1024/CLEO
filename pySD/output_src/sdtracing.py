'''
----- CLEO -----
File: sdtracing.py
Project: output_src
Created Date: Tuesday 24th October 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
Last Modified: Tuesday 24th October 2023
Modified By: CB
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
Copyright (c) 2023 MPI-M, Clara Bayley
-----
File Description:
functions to extract attribute data for
specifc superdroplets based on the sdIds
e.g. for tracing their trajectories
'''

import numpy as np
import awkward as ak

def superdrop_attribute(sddata, Id, attr):
  '''selects attribute from sddata belonging
  to superdroplet with identitiy 'Id'
  at every output time '''

  bools = ak.Array(sddata.sdId==Id) # True/False id is found in sdId at each time
  attr4Id = sddata[attr][bools] # attribute where sdId = Id

  num = ak.num(attr4Id) # at each time, number of positions where sdId is found (should be 0 or 1)
  if any(num[num!=1]):
    errmsg = "attribute timeseries has times when more"+\
      " than one sdId==Id. num should be list of either 1 or 0"+\
      " (for Id found in sddata at given time or not)"
    raise ValueError(errmsg)

  attr4Id = ak.where(num==0, ak.Array([[np.nan]]), attr4Id) # replace empty values with np.nan

  return ak.flatten(attr4Id, axis=1) # remove excess dimension

def attrtimeseries_for_superdropssample(sddata, attr,
                                        ndrops2sample=0,
                                        minid=0, maxid=0,
                                        ids=[]):
  ''' returns 2D array with dimensions [time, SD]
  containing attribute data over time for 'ndrops'
  randomly selected from superdrops with id in
  range [minid, maxid] '''

  if ids == []:
    population = list(range(minid, maxid, 1))
    sampled_ids = random.sample(population, ndrops2sample)
  else:
    sampled_ids = ids

  ndrops_attr = []
  for id in sampled_ids: 
    attr_timeseries = attrtimeseries_for_1superdrop(sddata, id, attr)
    ndrops_attr.append(attr_timeseries)
  
  return np.asarray(ndrops_attr).T

def attr_at_times(attrdata, time, times2sel):
  '''selects attribute at given times
   (for all superdroplets in sddata)'''

  selected_attr = [] # list containing an attribute at selected times
  for t in times2sel:
    ind = np.argmin(abs(time-t))
    selected_attr.append(attrdata[ind]) 
  
  return ak.Array(selected_attr)

def attrs_at_times(sddata, time, times2sel, attrs2sel):
  '''selects attributes at given times from
  sddata (for all superdroplets in sddata)'''

  selected_data = {} # dict containting selected attributes at selected times
  
  for attr in attrs2sel:
    
    selattr_data = attr_at_times(sddata[attr], time, times2sel)
    selected_data[attr] = selattr_data
  
  return selected_data

def surfaceprecip_fromdataset(dataset, gbxs):

    ds = get_rawdataset(dataset)

    sdindex = ak.unflatten(ds["sdindex"].values, ds["raggedcount"].values)
    radius = ak.unflatten(ds["radius"].values, ds["raggedcount"].values)
    eps = ak.unflatten(ds["eps"].values, ds["raggedcount"].values)

    r3sum = []
    for ti in range(ds.time.shape[0]-1):
        sd_ti, r_ti, eps_ti = sdindex[ti], radius[ti], eps[ti]
        sds_gone = set(sd_ti) - set(sdindex[ti+1]) # set of SD indexes that have left domain during timestep ti -> ti+1
        isgone = np.where(np.isin(sd_ti, list(sds_gone))) # indexes in ragged arrays of SDs that leave during timestep ti -> ti+1
        r3sum.append(np.dot(r_ti[isgone]**3, eps_ti[isgone])) # sum of (real) droplet radii^3 that left domain [microns^3]
    precipvol = 4/3 * np.pi * np.asarray(r3sum) / (1e18) # volume of water that left domain [m^3]

    domainy = np.amax(gbxs["yhalf"]) - np.amin(gbxs["yhalf"]) # [m]
    domainx = np.amax(gbxs["xhalf"]) - np.amin(gbxs["xhalf"]) # [m]
    deltat = np.diff(ds["time"].values) / 60 / 60 # [hrs]
    preciprate = precipvol * 1000 / (domainx * domainy) / deltat # [mm/hr]

    precipaccum = np.cumsum(preciprate * deltat) # [mm]
    preciprate = np.insert(preciprate, 0, 0) # at t=0, precip rate = 0
    precipaccum = np.insert(precipaccum, 0, 0) # at t=0, accumulated precip = 0

    return preciprate, precipaccum # [mm/hr] , [mm]