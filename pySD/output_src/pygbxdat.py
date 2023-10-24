'''
----- CLEO -----
File: pygbxdat.py
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
for reading data from gridbox (gbx)
boudnaries binary inputs files
'''

import numpy as np
from ..gbxboundariesbinary_src import read_gbxboundaries as rgrid

def get_gridboxes(gridfile, COORD0, isprint=True):

  gbxbounds, ndims =  rgrid.read_dimless_gbxboundaries_binary(gridfile,
                                                              COORD0=COORD0,
                                                              return_ndims=True,
                                                              isprint=isprint) 
  zhalf, xhalf, yhalf = rgrid.halfcoords_from_gbxbounds(gbxbounds,
                                                        isprint=isprint)
  domainvol, gbxvols, ngrid = rgrid.domaininfo(gbxbounds, isprint=isprint)
 
  gbxs = {
    "ngrid": ngrid, # number of gridboxes 
    "ndims": np.flip(ndims), # dimensions (no. gridboxes in [y,x,z] direction)
    "domainvol": domainvol,
    "gbxvols": gbxvols, # list of vols of each gbx 
    
    "zhalf": zhalf, # half cell coords (boundaries)
    "zfull": rgrid.fullcell(zhalf), # full cell coords (centres)
    
    "xhalf": xhalf, # half cell coords (boundaries)
    "xfull": rgrid.fullcell(xhalf), # full cell coords (centres)
    
    "yhalf": yhalf, # half cell coords (boundaries)
    "yfull": rgrid.fullcell(yhalf), # full cell coords (centres)
  }
  gbxs["gbxvols"] = np.reshape(gbxs["gbxvols"], gbxs["ndims"])

  xxh, zzh = np.meshgrid(gbxs["xhalf"], gbxs["zhalf"], indexing="ij") # dims [xdims, zdims] [m]
  xxf, zzf = np.meshgrid(gbxs["xfull"], gbxs["zfull"], indexing="ij") # dims [xdims, zdims] [m]
  gbxs["xxh"], gbxs["zzh"] = xxh, zzh
  gbxs["xxf"], gbxs["zzf"] = xxf, zzf
  
  return gbxs 