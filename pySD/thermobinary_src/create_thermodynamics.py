import numpy as np
from os.path import isfile
from .. import cxx2py, writebinary
from ..gbxboundariesbinary_src.read_gbxboundaries import read_dimless_gbxboundaries_binary

def get_Mrratio_from_constsfile(constsfile):
  
  consts = cxx2py.read_cpp_into_floats(constsfile, False)[0]
  moreconsts = cxx2py.derive_more_floats(consts, False)

  return moreconsts["Mr_ratio"]

def saturation_press(TEMP):
  ''' Calculate the equilibrium vapor pressure of water over
  liquid water ie. the saturation pressure (psat) given the
  temperature [K]. Equation from Bjorn Steven's "make_tetens"
  function in module "moist_thermodynamics.saturation_vapour_pressures"
  available on gitlab. Original paper "Murray, F. W. On the
  Computation of Saturation Vapor Pressure. Journal of Applied
  Meteorology and Climatology 6, 203–204 (1967). '''
  
  Aconst = 17.4146
  Bconst = 33.639
  TREF = 273.16 # Triple point temperature [K] of water
  PREF = 611.655 # Triple point pressure [Pa] of water

  if (TEMP <= 0.0):
    raise ValueError("psat ERROR: T must be larger than 0K."+\
                     " T = " + str(TEMP))

  return PREF * np.exp(Aconst * (TEMP - TREF) / (TEMP - Bconst)) # [Pa]


def relh2qvap(press, temp, relh, Mr_ratio):
  ''' convert relative humidity [%] (relh) into vapour mass
  mixing ratio (qvap) given ambient temperature and pressure
  and ratio of molecular masses: vapour/air '''
  
  vapourpress = saturation_press(temp) * relh / 100.0  # [Pa]
  
  qvap = Mr_ratio * vapourpress / (press - vapourpress) # dimensionless [Kg/kg]

  return qvap

  
class ConstUniformThermo:
  ''' create thermodyanmcis thats constant in time 
  and uniform throughout the domain '''

  def __init__(self, PRESS, TEMP, relh,
              qcond, WVEL, UVEL, VVEL,
              constsfile):
    self.PRESS = PRESS                        # pressure [Pa]
    self.TEMP = TEMP                          # temperature [T]
    
    Mr_ratio = get_Mrratio_from_constsfile(constsfile)
    self.qvap = relh2qvap(PRESS, TEMP, 
                          relh, Mr_ratio)     # water vapour content []
    
    self.qcond = qcond                        # liquid water content []
    self.WVEL = WVEL                          # vertical (z) velocity [m/s]
    self.UVEL = UVEL                          # horizontal x velocity [m/s]
    self.VVEL = VVEL                          # horizontal y velocity [m/s]

  def generate_thermo(self, ngrid, ntime):

    THERMODATA = {
      "PRESS": np.full(ngrid*ntime, self.PRESS),
      "TEMP": np.full(ngrid*ntime, self.TEMP),
      "qvap": np.full(ngrid*ntime, self.qvap),
      "qcond": np.full(ngrid*ntime, self.qcond), 
      "WVEL": np.full(ngrid*ntime, self.WVEL),
      "UVEL": np.full(ngrid*ntime, self.UVEL),
      "VVEL": np.full(ngrid*ntime, self.VVEL) 
    }

    return THERMODATA

def thermoinputsdict(configfile, constsfile):
  ''' create values from constants file & config file
  required as inputs to create initial 
  superdroplet conditions '''

  consts = cxx2py.read_cpp_into_floats(constsfile)[0]
  moreconsts = cxx2py.derive_more_floats(consts)
  config = cxx2py.read_configtxt_into_floats(configfile)[0]

  inputs = {
    # for creating thermodynamic profiles
    "G": consts["G"],
    "CP_DRY": consts["CP_DRY"],
    "RHO_DRY": consts["RHO_DRY"],               # dry air density [Kg/m^3]
    "RGAS_DRY": moreconsts["RGAS_DRY"],
    "Mr_ratio": moreconsts["Mr_ratio"],
    "COUPLTSTEP": config["COUPLTSTEP"],
    "T_END": config["T_END"],

    # for de-dimensionalising attributes
    "W0": consts["W0"],
    "P0": consts["P0"],
    "TEMP0": consts["TEMP0"],
    "RHO0": moreconsts["RHO0"],               # characteristic density scale [Kg/m^3]
    "CP0": moreconsts["CP0"],
    "COORD0": moreconsts["COORD0"]            # z coordinate lengthscale [m]
  }

  return inputs

def dimless_thermodynamics(THERMO, inputs):

  sfs = [inputs["P0"], inputs["TEMP0"],
                   1.0, 1.0] + [inputs["W0"]]*3 # scale_factors to de-dimensionalise data

  thermodata = {
      "press": THERMO["PRESS"] / sfs[0],
      "temp":THERMO["TEMP"] / sfs[1],
      "qvap":THERMO["qvap"] / sfs[2],
      "qcond": THERMO["qcond"] / sfs[3],
      "wvel":THERMO["WVEL"] / sfs[4],
      "uvel":THERMO["UVEL"] / sfs[5],
      "vvel":THERMO["VVEL"] / sfs[6]
    }
  
  return thermodata, sfs
  

def write_thermodynamics_binary(thermofile, thermogen, configfile,
                                constsfile, gridfile):
  
  if not isfile(gridfile):
    errmsg = "gridfile not found, but must be"+\
              " created before initSDsfile can be"
    raise ValueError(errmsg)

  inputs = thermoinputsdict(configfile, constsfile)
  gbxbounds = read_dimless_gbxboundaries_binary(gridfile,
                                                COORD0=inputs["COORD0"])
  
  ngridboxes = len(gbxbounds.keys())
  ntime = round(inputs["T_END"]/inputs["COUPLTSTEP"])+1
  thermodata = thermogen.generate_thermo(ngridboxes, ntime)

  print(thermodata["PRESS"], thermodata["TEMP"], thermodata["qvap"])
  print(saturation_press(thermodata["TEMP"][0]))

  thermodata, scale_factors = dimless_thermodynamics(thermodata, inputs)
