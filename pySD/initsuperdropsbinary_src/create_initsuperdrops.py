import numpy as np
from os.path import isfile
from .. import cxx2py, writebinary
from ..gbxboundariesbinary_src.read_gbxboundaries import read_dimless_gbxboundaries_binary

class ManyInitAttrs:
    '''store for lists of each attribute for all superdroplets ''' 
    def __init__(self):
        self.sd_gbxindex = []
        self.eps = []
        self.radius = []
        self.m_sol = []
        self.coord3 = []
        self.coord1 = []
        self.coord2 = []
    
    def set_attrlists(self, a, b, c, 
                      d, e, f, g):
        self.sd_gbxindex = a
        self.eps = b
        self.radius = c
        self.m_sol = d
        self.coord3 = e 
        self.coord1 = f 
        self.coord2 = g 
        
    def extend_attrlists_fromlists(self, a, b, c, 
                                   d, e, f, g):
        self.sd_gbxindex.extend(a)
        self.eps.extend(b)
        self.radius.extend(c)
        self.m_sol.extend(d)
        self.coord3.extend(e)
        self.coord1.extend(f) 
        self.coord2.extend(g) 
   
    def extend_attrlists(self, mia):
        ''' use an instance of ManyInitAttrs (mia) to
        extend lists in this instance '''
        self.sd_gbxindex.extend(mia.sd_gbxindex)
        self.eps.extend(mia.eps)
        self.radius.extend(mia.radius)
        self.m_sol.extend(mia.m_sol)
        self.coord3.extend(mia.coord3)
        self.coord1.extend(mia.coord1)
        self.coord2.extend(mia.coord2)   

def initSDsinputsdict(configfile, constsfile):
  ''' create values from constants file & config file
  required as inputs to create initial 
  superdroplet conditions '''

  consts = cxx2py.read_cpp_into_floats(constsfile)[0]
  moreconsts = cxx2py.derive_more_floats(consts)
  config = cxx2py.read_configtxt_into_floats(configfile)[0]

  inputs = {
    # for creating SD attribute distirbutions
    "SDnspace": int(config["SDnspace"]),
    "RHO_SOL": consts["RHO_SOL"],               # solute density [Kg/m^3]

    # for de-dimensionalising attributes
    "R0":consts["R0"],                          # droplet radius lengthscale [m]
    "RHO0": moreconsts["RHO0"],                 # characteristic density scale [Kg/m^3]
    "COORD0": moreconsts["COORD0"],                              # z coordinate lengthscale [m]
  }

  inputs["MASS0"] = inputs["RHO0"] * (inputs["R0"]**3)

  return inputs

def is_sd_gbxindex_correct(gridboxbounds, coord3,
                           sd_gbxindex, gbxindex, COORD0):

  coord3 = coord3*COORD0 #[m]
  
  errmsg = None
  if (coord3 < gridboxbounds[0]).any():
    errmsg = "superdroplet coord3 below lower"+\
              " limit of gridbox it's associated with"
  elif (coord3 >= gridboxbounds[1]).any(): 
    errmsg = "superdroplet coord3 above upper"+\
              " limit of gridbox it's associated with" 
  elif (sd_gbxindex != gbxindex).any():
    errmsg = "superdroplet gridbox index not same as"+\
              " gridbox it should be associated with"
  
  if errmsg:
    raise ValueError(errmsg)

def dimless_superdropsattrs(nsupers, initattrsgen, inputs, gbxindex,
                            gridboxbounds, NUMCONC):
    ''' use superdroplet attribute generator "initattrsgen"
    and settings from config and consts files to
    make dimensionless superdroplet attributes'''
    
    # generate attributes
    sd_gbxindex = [gbxindex]*nsupers
    eps, radius, m_sol = initattrsgen.generate_attributes(nsupers, 
                                                          inputs["RHO_SOL"],
                                                          NUMCONC,
                                                          gridboxbounds) 
    coord3, coord1, coord2 = initattrsgen.generate_coords(nsupers,
                                                          inputs["SDnspace"],
                                                          gridboxbounds)

    # de-dimsionalise attributes
    radius = radius / inputs["R0"]
    m_sol = m_sol / inputs["MASS0"]
    coord3 = coord3 / inputs["COORD0"]
    coord1 = coord1 / inputs["COORD0"]
    coord2 = coord2 / inputs["COORD0"]

    attrs4gbx = ManyInitAttrs() 
    attrs4gbx.set_attrlists(sd_gbxindex, eps, radius, m_sol,
                            coord3, coord1, coord2)

    is_sd_gbxindex_correct(gridboxbounds, coord3, sd_gbxindex,
                           gbxindex, inputs["COORD0"])

    return attrs4gbx

def create_allsuperdropattrs(nsupersdict, initattrsgen,
                             gbxbounds, inputs, NUMCONC):
  ''' returns lists for attributes of all SDs in domain called attrs'''

  attrs = ManyInitAttrs() # lists of attrs for SDs in domain

  for gbxindex, gridboxbounds in gbxbounds.items():

    nsupers = nsupersdict[gbxindex]
    attrs4gbx = dimless_superdropsattrs(nsupers, initattrsgen, inputs,
                                        gbxindex, gridboxbounds, NUMCONC) # lists of attrs for SDs in gridbox

    attrs.extend_attrlists(attrs4gbx)

  return attrs

def set_arraydtype(arr, dtype):
   
  og = type(arr[0])
  if og != dtype: 
    arr = np.array(arr, dtype=dtype)

    warning = "WARNING! dtype of attributes is being changed!"+\
                " from "+str(og)+" to "+str(dtype)
    raise ValueError(warning) 

  return arr

def ctype_compatible_attrs(attrs):
  ''' make list from arrays of SD attributes that are compatible
  with c type expected by SDM e.g. unsigned long ints for eps,
  doubles for radius and m_sol'''   

  datatypes = [np.uintc, np.uint, np.double, np.double]
  datatypes += [np.double]*3 # coords datatype
  
  attrs.sd_gbxindex = list(set_arraydtype(attrs.sd_gbxindex, datatypes[0]))
  attrs.eps = list(set_arraydtype(attrs.eps, datatypes[1]))
  attrs.radius = list(set_arraydtype(attrs.radius, datatypes[2]))
  attrs.m_sol = list(set_arraydtype(attrs.m_sol, datatypes[3]))
  
  datalist = attrs.sd_gbxindex + attrs.eps + attrs.radius + attrs.m_sol
  
  if any(attrs.coord3):
    # make coord3 compatible if there is data for it (>= 1-D model)
    attrs.coord3 = list(set_arraydtype(attrs.coord3, datatypes[4]))
    datalist += attrs.coord3

    if any(attrs.coord1):
      # make coord1 compatible if there is data for it and coord3 (>= 2-D model)
      attrs.coord1 = list(set_arraydtype(attrs.coord1, datatypes[4]))
      datalist += attrs.coord1
    
      if any(attrs.coord2):
        # make coord2 compatible if there is data for it and coord3 and coord2 (>= 3-D model)
        attrs.coord2 = list(set_arraydtype(attrs.coord2, datatypes[4]))
        datalist += attrs.coo2d3

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

def nsupers_pergridboxdict(nsupers, gbxbounds):
  
  if type(nsupers) == int:
    nsupersdict = {}
    for key in gbxbounds.keys():
      nsupersdict[key] = nsupers
    return nsupersdict

  elif type(nsupers) == dict:
    print(nsupers.keys())
    if nsupers.keys() != gbxbounds.keys():
      errmsg = "keys for nsupers dict don't match gridbox indexes"
      raise ValueError(errmsg)
    else:
      return nsupers
      nsupersdict == nsupers

  else:
    errmsg = "nsupers must be either dict or int"
    raise ValueError(errmsg)


def write_initsuperdrops_binary(initSDsfile, initattrsgen, configfile, 
                                constsfile, gridfile, nsupers, NUMCONC):
  ''' de-dimensionalise attributes in initattrsgen and then write to 
  to a binary file, "initSDsfile", with some metadata '''
  
  if not isfile(gridfile):
    errmsg = "gridfile not found, but must be"+\
              " created before initSDsfile can be"
    raise ValueError(errmsg)

  inputs = initSDsinputsdict(configfile, constsfile)
  gbxbounds = read_dimless_gbxboundaries_binary(gridfile,
                                                COORD0=inputs["COORD0"])
  nsupersdict = nsupers_pergridboxdict(nsupers, gbxbounds) 
  
  attrs = create_allsuperdropattrs(nsupersdict, initattrsgen,
                                   gbxbounds, inputs, NUMCONC) 
  
  ndata = [len(dt) for dt in [attrs.sd_gbxindex, attrs.eps,
                              attrs.radius, attrs.m_sol, attrs.coord3,
                              attrs.coord1, attrs.coord2]]
  data, datatypes = ctype_compatible_attrs(attrs) 
  check_datashape(data, ndata)

  units = [b' ', b' ', b'm', b'g']
  units += [b'm']*3 # coords units
  
  scale_factors = [1.0, 1.0, inputs["R0"], inputs["MASS0"]]
  scale_factors += [inputs["COORD0"]]*3 # coords scale factors
  scale_factors = np.asarray(scale_factors, dtype=np.double)

  metastr = 'Variables in this file are Superdroplet attributes:'
  if initattrsgen.coord3gen: 
    if initattrsgen.coord1gen:
      if initattrsgen.coord2gen: 
        metastr += ' [sd_gbxindex, eps, radius, m_sol, coord3, coord1, coord2]'
      else:
        metastr += ' [sd_gbxindex, eps, radius, m_sol, coord3, coord1]'
    else:
      metastr += ' [sd_gbxindex, eps, radius, m_sol, coord3]'
  else:
    metastr += ' [sd_gbxindex, eps, radius, m_sol]'
  
  writebinary.writebinary(initSDsfile, data, ndata, datatypes,
                          units, scale_factors, metastr)