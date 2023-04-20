import numpy as np
from .. import cxx2py

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
      "WVEL": np.array([]), 
      "UVEL": np.array([]),
      "VVEL": np.array([])
    }

    if self.WVEL != None:
      THERMODATA["WVEL"] =  np.full(ngrid*ntime, self.WVEL)

      if self.UVEL != None:
        THERMODATA["UVEL"] =  np.full(ngrid*ntime, self.UVEL)

        if self.VVEL != None:
          THERMODATA["VVEL"] =  np.full(ngrid*ntime, self.VVEL)

    # THERMODATA["PRESS"][0:ngrid] = 800 # makes all gbxs a t=0 have P=800Pa
    # THERMODATA["PRESS"][::ngrid] = 800 # makes 0th gbx at all times have P=800Pa

    return THERMODATA

