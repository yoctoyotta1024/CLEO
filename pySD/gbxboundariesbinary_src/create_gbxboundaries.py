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

def linearlyspaced_halfcoords(min, max, delta):
  ''' returns linearly spaced half coords ie. 1D 
  gridbox boundaries given the maximum, minimum coord
  and the spacing '''

  return np.arange(min, max+delta, delta, dtype=np.double)

def linearlyspaced_zhalf(zmin, zmax, zdelta):
  ''' returns linearly spaced z half coordinates
  given the z coord limits and spacing. returned
  zhalf goes from largest to smallest z coord '''

  zhalf = np.flip(linearlyspaced_halfcoords(zmin, zmax, zdelta))

  return zhalf

def halfcoords_from_coordlims(grid, coord):
  ''' return half coordinates (ie. gridbox boundaries)
  given 'grid' list of [min coord, max coord, delta coord]
  for x, y or z coord '''

  if (coord not in ["z", "x", "y"]):
    raise ValueError("coord should be x, y or z")

  elif coord == "z":
    return linearlyspaced_zhalf(grid[0], grid[1], grid[2]) 
  
  elif coord == "x" or coord == "y":
    return linearlyspaced_halfcoords(grid[0], grid[1], grid[2]) 

def check_halfcoords(grid, coord):
  ''' return half coordinates (ie. gridbox boundaries)
  given array containing 
  for x, y or z coord '''
  
  if (coord not in ["z", "x", "y"]):
    raise ValueError("coord should be x, y or z")

  elif coord == "z":
    # check z grid is at least 1 cell with
    # strictly monotonically decreasing boundaries
    dz = np.diff(grid)
    criteria = (len(grid) >= 2 and np.all(dz < 0))
    
  elif coord == "x" or coord == "y":
    # check that x or y grid is single cell with
    # strictly monotonically increasing bounds
    criteria = (len(grid) == 2 and np.all(np.diff(grid) > 0))

  if criteria != True:
    raise ValueError(str(grid)+" does not meet criteria for "+coord+" halfcoords")

def get_dimless_halfcoords(grid, coord, COORD0):
  ''' given a list or array, return dimensionless halfcoords.
  For list call "halfcoords_from_coordlims" before
  de-dimensionalising. For array simply return array / COORD0'''

  if (type(grid) not in [list, np.ndarray]):
    raise ValueError("input "+coord+" grid is neither list or array ") 

  elif type(grid) == list:
    grid = halfcoords_from_coordlims(grid, coord)

  elif type(grid) == np.ndarray:
    grid = np.array(grid, dtype=np.double)
   
  return grid / COORD0

def dimless_gridboxboundaries(zgrid, xgrid, ygrid, COORD0):
  ''' use zgrid, xgrid and ygrid lists or arrays to create half coords
  of gridboxes (ie. their boundaries). Return single list of zhalf, 
  then xhalf, then yhalf coords and a list of len(3) which states how
  many of each z, x and y coords there are'''

  allhalfcoords = [] # z, x and y half coords in 1 continuous list
  num_halfcoords = [] # number of z, x and y half coords
  for grid, coord in zip([zgrid, xgrid, ygrid], ["z", "x", "y"]):
    
    halfcoords = get_dimless_halfcoords(grid, coord, COORD0)
    
    check_halfcoords(halfcoords, coord)
    
    allhalfcoords.extend(halfcoords)
    num_halfcoords.append(len(halfcoords))

  return allhalfcoords, num_halfcoords

def ctype_compatible_gridboxboundaries(data):
  ''' check type of gridbox boundaries data is compatible
  with c type double. If not, change type and raise error '''

  og = type(data[0])
  if og != np.double: 
    
    data = [np.double(d) for d in data]
    
    warning = "WARNING! dtype of gbx boundaries is being changed!"+\
              " from "+str(og)+" to np.double"
    raise ValueError(warning)
    #print(warning)
  
  datatypes = [np.double]*3

  return data, datatypes

def write_gridboxboundaries_binary(filename, zgrid, xgrid, ygrid, constsfile):
  ''' zgrid, xgrid and ygrid can either be list of
  [min coord, max coord, delta] or they can be arrays of
  their half coordinates (ie. gridbox boundaries). If the former,
  first create half coordinates. In both cases, de-dimensionalise
  and write the boundaries to a binary file, "filename" '''

  COORD0 = get_COORD0_from_constsfile(constsfile)

  halfcoordsdata, nhalfcoords = dimless_gridboxboundaries(zgrid, xgrid, 
                                                            ygrid, COORD0) 

  metastr = 'Variables in this file are coordinates of'+\
            ' [zhalf, xhalf, yhalf] gridbox boundaries'
  
  halfcoordsdata, datatypes = ctype_compatible_gridboxboundaries(halfcoordsdata) 
  scale_factors = np.array([COORD0, COORD0, COORD0], dtype=np.double) 
  units = [b"m"]*3
  writebinary.writebinary(filename, halfcoordsdata, nhalfcoords, 
                datatypes, units, scale_factors, metastr)