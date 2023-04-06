import numpy as np

from .. import cxx2py, writebinary

def get_COORD0_from_constsfile(constsfile, returnconsts=False):
  ''' create values from constants file & config file
  required as inputs to create initial 
  superdroplet conditions '''

  consts = cxx2py.read_cpp_into_floats(constsfile)[0]
  COORD0 = consts["TIME0"]*consts["W0"]
  
  if returnconsts:
    return COORD0, consts
  else:
    return COORD0

def linearlyspaced_halfcoords(grid):
  ''' returns linearly spaced half coords ie. 1D 
  gridbox boundaries given 'grid' list of
  [min coord, max coord, delta coord] for x, y or z'''

  min, max, delta = grid

  return np.arange(min, max+delta, delta, dtype=np.double)

def get_dimless_halfcoords(grid, coord, COORD0):
  ''' given a list or array, return dimensionless halfcoords.
  For list call "halfcoords_from_coordlims" before
  de-dimensionalising. For array simply return array / COORD0'''

  if (type(grid) not in [list, np.ndarray]):
    raise ValueError("input "+coord+" grid is neither list or array ") 

  elif type(grid) == list:
    grid = linearlyspaced_halfcoords(grid) 

  elif type(grid) == np.ndarray:
    grid = np.array(grid, dtype=np.double)
   
  return grid / COORD0

def check_halfcoords(grid, coord):
  ''' check that x , y or z grid limits are for at least
  1 cell with strictly monotonically increasing bounds '''
  
  if (coord not in ["z", "x", "y"]):
    raise ValueError("coord should be x, y or z")
  
  criteria = (len(grid) >= 2 and np.all(np.diff(grid) > 0))
  if criteria != True:
    errmsg = str(grid)+" does not meet criteria for "+coord+" halfcoords"
    raise ValueError(errmsg)

def gridboxboundaries_from_halfcoords(allhalfcoords):
  ''' returns gbxbounds dictionary. Each key of gbxbounds is a gridbox's
  index and the corresponding value is the gridbox's
  [zmin, zmax, xmin, xmax, ymin, ymax] boundaries'''

  zhalfs = allhalfcoords[0]
  xhalfs = allhalfcoords[1]
  yhalfs = allhalfcoords[2]

  gbxbounds, ii = {}, 0
  for j in range(len(yhalfs)-1):
    for i in range(len(xhalfs)-1):
      for k in range(len(zhalfs)-1):
        zbounds = [zhalfs[k], zhalfs[k+1]]
        xbounds = [xhalfs[i], xhalfs[i+1]]
        ybounds = [yhalfs[j], yhalfs[j+1]]
        gbxbounds[ii] = zbounds + xbounds + ybounds
        ii+=1
 
  ngridboxes = len(gbxbounds)
  print("created boundaries for",ngridboxes,"gridboxes")

  return gbxbounds, ngridboxes

def dimless_gridboxboundaries(zgrid, xgrid, ygrid, COORD0):
  ''' use zgrid, xgrid and ygrid lists or arrays to create half coords
  of gridboxes (ie. their boundaries). Return single list of zhalf, 
  then xhalf, then yhalf coords and a list of len(3) which states how
  many of each z, x and y coords there are'''

  allhalfcoords = [] # z, x and y half coords in 1 continuous list
  for grid, coord in zip([zgrid, xgrid, ygrid], ["z", "x", "y"]):
    
    halfcoords = get_dimless_halfcoords(grid, coord, COORD0)
    check_halfcoords(halfcoords, coord)
  
    allhalfcoords.append(halfcoords)
  
  gbxbounds, ngridboxes = gridboxboundaries_from_halfcoords(allhalfcoords)
  
  gbxindicies = np.array(list(gbxbounds.keys()), dtype=np.uintc).flatten()
  gbxboundsdata = np.array(list(gbxbounds.values()), dtype=np.double).flatten()
  
  return gbxindicies, gbxboundsdata, ngridboxes

def set_arraydtype(arr, dtype):
   
  og = type(arr[0])
  if og != dtype: 
    arr = np.array(arr, dtype=dtype)

    warning = "WARNING! dtype of attributes is being changed!"+\
                " from "+str(og)+" to "+str(dtype)
    raise ValueError(warning) 

  return arr

def ctype_compatible_gridboxboundaries(idxs, bounds):
  ''' check type of gridbox boundaries data is compatible
  with c type double. If not, change type and raise error '''

  datatypes = [np.uintc, np.double]

  idxs = list(set_arraydtype(idxs, datatypes[0]))
  bounds = list(set_arraydtype(bounds, datatypes[1]))

  datalist = idxs + bounds

  return datalist, datatypes

def write_gridboxboundaries_binary(gridfile, zgrid, xgrid, ygrid, constsfile):
  ''' zgrid, xgrid and ygrid can either be list of
  [min coord, max coord, delta] or they can be arrays of
  their half coordinates (ie. gridbox boundaries). If the former,
  first create half coordinates. In both cases, de-dimensionalise
  and write the boundaries to a binary file, "filename" '''

  COORD0 = get_COORD0_from_constsfile(constsfile)

  gbxindicies, gbxboundsdata, ngridboxes = dimless_gridboxboundaries(zgrid, xgrid,
                                                                     ygrid, COORD0) 
  
  metastr = 'Variables in this file are '+str(ngridboxes)+' gridbox'+\
            ' indicies followed by the [zmin, zmax, xmin, xmax, ymin, ymax]'+\
            ' coordinates for each gridbox\'s boundaries'
  
  ndata = [len(dt) for dt in [gbxindicies, gbxboundsdata]]
  data, datatypes = ctype_compatible_gridboxboundaries(gbxindicies, gbxboundsdata) 
  scale_factors = np.array([COORD0, COORD0, COORD0], dtype=np.double) 
  units = [b' ', b"m"]

  writebinary.writebinary(gridfile, data, ndata, datatypes,
                          units, scale_factors, metastr)