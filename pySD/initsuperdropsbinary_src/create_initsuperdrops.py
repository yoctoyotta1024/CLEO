import numpy as np
from .. import cxx2py, writebinary

def initSDsinputsdict(configfile, constsfile):
  ''' create values from constants file & config file
  required as inputs to create initial 
  superdroplet conditions '''

  consts = cxx2py.read_cpp_into_floats(constsfile)[0]
  moreconsts = cxx2py.derive_more_floats(consts)
  config = cxx2py.read_configtxt_into_floats(configfile)[0]

  inputs = {
    
    # for creating SD attribute distirbutions
    "nsupers": int(config["NSUPERS"]),          # no. of distinct superdrops (different initial radii (evenly spaced between ln(rspan))
    "SDnspace": int(config["SDnspace"]),
    "RHO_SOL": consts["RHO_SOL"],               # solute density [Kg/m^3]

    # for de-dimensionalising attributes
    "R0":consts["R0"],                          # droplet radius lengthscale [m]
    "RHO0": moreconsts["RHO0"],                 # characteristic density scale [Kg/m^3]
    "COORD0": moreconsts["COORD0"],                              # z coordinate lengthscale [m]
  }

  inputs["MASS0"] = inputs["RHO0"] * (inputs["R0"]**3)

  return inputs
    

def dimless_superdropsattributes(initattrs, inputs):
    ''' use superdroplet attribute generator "initattrs"
    and settings from config and consts files to
    make dimensionless superdroplet attributes'''
    
    # generate attributes
    epss, radii, m_sols, coord3s = initattrs.generate_attributes(inputs["nsupers"], 
                                                              inputs["RHO_SOL"],
                                                              inputs["SDnspace"])
    
    # de-dimsionalise attributes
    radii = radii / inputs["R0"]
    m_sols = m_sols / inputs["MASS0"]
    coord3s = coord3s / inputs["COORD0"]

    return epss, radii, m_sols, coord3s

def set_arraydtype(arr, dtype):
   
  og = type(arr[0])
  if og != dtype: 
    arr = np.array(arr, dtype=dtype)

    warning = "WARNING! dtype of attributes is being changed!"+\
                " from "+str(og)+" to "+str(dtype)
    raise ValueError(warning) 

  return arr

def ctype_compatible_attrs(epss, radii, m_sols, coord3s):
  ''' make list from arrays of SD attributes that are compatible
  with c type expected by SDM e.g. unsigned long ints for eps,
  doubles for radius and m_sol'''   

  datatypes = [np.uint, np.double, np.double, np.double]
  
  epss = list(set_arraydtype(epss, datatypes[0]))
  radii = list(set_arraydtype(radii, datatypes[1]))
  m_sols = list(set_arraydtype(m_sols, datatypes[2]))
  
  datalist = epss + radii + m_sols
  
  if any(coord3s):
    coord3s = list(set_arraydtype(coord3s, datatypes[3]))
    datalist += coord3s

  return datalist, datatypes

def check_datashape(data, ndata):
  ''' make sure each superdroplet attribute in data has length stated
  in ndata and that this length is compatible with the nummber of
  attributes and superdroplets expected given ndata'''
  
  if any([n != len(data) / len(ndata) for n in ndata]):
    
    print("\n------ WARNING! -----\n",
          "not all variables in data are same length, ndata = ",
          ndata, "\n---------------------\n")
    
    if len(data) != np.sum(ndata): 
      err = "inconsistent dimensions of data: "+str(np.shape(data))+", and"+\
            " data per attribute: "+str(ndata)+". data should be 1D with"+\
            " shape: num_attributes * nsupers. nata should be list of"+\
            " [nsupers]*num_attributes."     
      raise ValueError(err)
                                                                                                                  
def write_initsuperdrops_binary(initSDsfile, initattrs, configfile, constsfile):

  ''' de-dimensionalise attributes in initattrs and then write to 
  to a binary file, "initSDsfile", with some metadata '''
  
  inputs = initSDsinputsdict(configfile, constsfile)

  epss, radii, m_sols, coord3s = dimless_superdropsattributes(initattrs, 
                                                              inputs) 
  ndata = [len(dt) for dt in [epss, radii, m_sols, coord3s]]
  data, datatypes = ctype_compatible_attrs(epss, radii, m_sols, coord3s) 
  check_datashape(data, ndata)

  units = [b' ', b'm', b'g', b'm']
  scale_factors = np.array([1.0, inputs["R0"], inputs["MASS0"], 
                           inputs["COORD0"]], dtype=np.double)

  if initattrs.coord3gen: 
    metastr = 'Variables in this file are Superdroplet attributes:'+\
            ' [eps, radius, m_sol, coord3]'
  else:
    metastr = 'Variables in this file are Superdroplet attributes:'+\
            ' [eps, radius, m_sol]'
  
  writebinary.writebinary(initSDsfile, data, ndata, datatypes,
                units, scale_factors, metastr)