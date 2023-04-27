import numpy as np
from .. import cxx2py
from ..gbxboundariesbinary_src.read_gbxboundaries import fullcoords_forallgridboxes

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

  if (np.any(TEMP <= 0.0)):
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

def sratio2qvap(sratio, press, temp, Mr_ratio):

  psat = saturation_press(temp)
  
  qvap = Mr_ratio * sratio 
  qvap = qvap / (press/psat - 1)

  return qvap

def divfree_flowfield2D(WMAX, rhotilda, Zlength, Xlength, 
                        ZCOORDS, XCOORDS):

    ztilda = np.pi * ZCOORDS / Zlength 
    xtilda = 2* np.pi * XCOORDS / Xlength 
    VELfactor = WMAX / rhotilda
        
    WVEL = 2 * VELfactor * np.sin(ztilda) * np.sin(xtilda)
    
    UVEL = VELfactor *  Xlength /  Zlength
    UVEL =  UVEL * np.cos(ztilda) * np.cos(xtilda)

    return WVEL, UVEL
  
class ConstUniformThermo:
  ''' create thermodyanmics that's constant in time 
  and uniform throughout the domain '''

  def __init__(self, PRESS, TEMP, qvap,
              qcond, WVEL, UVEL, VVEL,
              relh=False, constsfile=""):
    self.PRESS = PRESS                        # pressure [Pa]
    self.TEMP = TEMP                          # temperature [T]
    
    if relh:
      Mr_ratio = get_Mrratio_from_constsfile(constsfile)
      self.qvap = relh2qvap(PRESS, TEMP, 
                            relh, Mr_ratio)     # water vapour content []
    else:
      self.qvap = qvap

    self.qcond = qcond                        # liquid water content []
    self.WVEL = WVEL                          # vertical (z) velocity [m/s]
    self.UVEL = UVEL                          # horizontal x velocity [m/s]
    self.VVEL = VVEL                          # horizontal y velocity [m/s]

  def generate_thermo(self, gbxbounds, ndims, ntime):

    ngridboxes = int(np.prod(ndims))

    THERMODATA = {
      "PRESS": np.full(ngridboxes*ntime, self.PRESS),
      "TEMP": np.full(ngridboxes*ntime, self.TEMP),
      "qvap": np.full(ngridboxes*ntime, self.qvap),
      "qcond": np.full(ngridboxes*ntime, self.qcond), 
      "WVEL": np.array([]), 
      "UVEL": np.array([]),
      "VVEL": np.array([])
    }

    if self.WVEL != None:
      THERMODATA["WVEL"] =  np.full(ngridboxes*ntime, self.WVEL)

      if self.UVEL != None:
        THERMODATA["UVEL"] =  np.full(ngridboxes*ntime, self.UVEL)

        if self.VVEL != None:
          THERMODATA["VVEL"] =  np.full(ngridboxes*ntime, self.VVEL)

    # THERMODATA["PRESS"][0:ngridboxes] = 800 # makes all gbxs a t=0 have P=800Pa
    # THERMODATA["PRESS"][::ngridboxes] = 800 # makes 0th gbx at all times have P=800Pa

    return THERMODATA

class SimpleThermo2Dflowfield:
  ''' create thermodyanmics that's constant in time 
  with (P,T,qc) uniform throughout the domain with relative humidity
  = 0.95 below zbase and a 2D (z,x) dependent flow field'''

  def __init__(self, PRESS, TEMP, qvapmethod, qcond, WMAX, Zlength,
               Xlength, VVEL, zbase, qparam, constsfile=''):
    
    self.PRESS = PRESS                        # pressure [Pa]
    self.TEMP = TEMP                          # temperature [T]
    self.qcond = qcond                        # liquid water content []
    
    # determine qvap above z (cloud) base
    self.zbase = zbase
    Mr_ratio = get_Mrratio_from_constsfile(constsfile)
    self.qvap_below = sratio2qvap(0.85, PRESS, TEMP,Mr_ratio) # supersat=0.85 below zbase
    if qvapmethod == "relh":
      self.qvap_above = relh2qvap(PRESS, TEMP, 
                            qparam, Mr_ratio)     # water vapour content []
    elif qvapmethod == "sratio":
      self.qvap_above = sratio2qvap(qparam, self.PRESS, self.TEMP, Mr_ratio)
    
    self.WMAX = WMAX  # max velocities constant
    self.Zlength = Zlength # wavelength of velocity modulation in z direction [m]
    self.Xlength = Xlength # wavelength of velocity modulation in x direction [m]
    self.VVEL = VVEL # horizontal (y) velocity

  def generate_qvap_profile(self, zfulls):

    qvap = np.where(zfulls >= self.zbase, self.qvap_above, self.qvap_below)

    return qvap

  def generate_thermo(self, gbxbounds, ndims, ntime):

      ngridboxes = int(np.prod(ndims))
      zfulls, xfulls, yfulls = fullcoords_forallgridboxes(gbxbounds, ndims)
      
      qvap = self.generate_qvap_profile(zfulls)

      THERMODATA = {
        "PRESS": np.full(ngridboxes*ntime, self.PRESS),
        "TEMP": np.full(ngridboxes*ntime, self.TEMP),
        "qvap": np.tile(qvap, ntime),
        "qcond": np.full(ngridboxes*ntime, self.qcond), 
        "WVEL": np.array([]), 
        "UVEL": np.array([]),
        "VVEL": np.array([])
      }

      if self.WMAX != None:
        WVEL, UVEL = divfree_flowfield2D(self.WMAX, 1.0, self.Zlength, 
                                          self.Xlength, zfulls, xfulls)
        THERMODATA["WVEL"] =  np.tile(WVEL, ntime)
        THERMODATA["UVEL"] =  np.tile(UVEL, ntime)

        if self.VVEL != None:
          THERMODATA["VVEL"] =  np.full(int(ngridboxes*ntime), self.VVEL)

      return THERMODATA
  
class ConstHydrostaticAdiabat:
  ''' create thermodyanmics that's constant in time 
  and in hydrostatic equillibrium with a dry adiabat 
  accounting for the mass of water vapour in the air.
  Equations derived from Arabas et al. 2015 (sect 2.1) '''

  def __init__(self, PRESS0, THETA, qvap, qcond, WMAX,
               Zlength, Xlength, VVEL, GRAVG, CP_DRY,
               RGAS_DRY, RGAS_V):
    
    ### parameters of profile ###
    self.PRESS0 = PRESS0 #pressure at z=0m [Pa]
    self.THETA = THETA # (constant) dry potential temperature [K]
    self.qvap = qvap # (constant) vapour mass mixing ratio []
    self.qcond = qcond # liquid mass mixing ratio []
    self.WMAX = WMAX  # max velocities constant
    self.Zlength = Zlength # wavelength of velocity modulation in z direction [m]
    self.Xlength = Xlength # wavelength of velocity modulation in x direction [m]
    self.VVEL = VVEL # horizontal (y) velocity

    ### constants ###
    self.GRAVG = GRAVG
    self.CP_DRY = CP_DRY
    self.RGAS_DRY = RGAS_DRY
    self.RGAS_V = RGAS_V
    self.RC_DRY = RGAS_DRY / CP_DRY
    self.RCONST = 1 + self.qvap * RGAS_V / RGAS_DRY
    self.P1000 = 100000 # P_1000 = 1000 hPa [Pa]
    
    alpha = PRESS0 / (self.RCONST * self.P1000) 
    self.TEMP0 = THETA * np.power(alpha, self.RC_DRY) #temperature at z=0m [K]

    beta = (1+self.qvap) / self.RCONST / self.RGAS_DRY
    self.RHO0 = beta * self.PRESS0 / self.TEMP0

  def hydrostatic_adiabatic_profile(self, ZCOORDS):
    ''' returns *profile* of density (not the density itself!)
    rho = rhoprofile^((1-RC_DRY)/RC_DRY) = profile^pow '''
    
    pow = 1/self.RC_DRY - 1
    
    Aa = (1+self.qvap)*np.power(self.P1000, self.RC_DRY)
    Aa = (self.THETA * self.RGAS_DRY / Aa) 
    Aconst = self.RCONST * np.power(Aa, (1 / (1-self.RC_DRY)))

    RHOconst = -1 * self.GRAVG * self.RC_DRY / Aconst
    RHOprofile = np.power(self.RHO0, 1/pow) + RHOconst * ZCOORDS # RHO0^pow

    return RHOprofile, Aconst

  def hydrostatic_adiabatic_thermo(self, ZCOORDS):
    
    RHOprof, Aconst = self.hydrostatic_adiabatic_profile(ZCOORDS)
    
    RHO = np.power(RHOprof, (1/self.RC_DRY - 1))

    PRESS = Aconst * np.power(RHOprof, 1/self.RC_DRY)
    
    TEMPconst = np.power(Aconst / (self.RCONST * self.P1000), self.RC_DRY)
    TEMP = self.THETA * TEMPconst * RHOprof
    
    return RHO, PRESS, TEMP

  def generate_thermo(self, gbxbounds, ndims, ntime):

    ngridboxes = int(np.prod(ndims))
    zfulls, xfulls, yfulls = fullcoords_forallgridboxes(gbxbounds, ndims)
    RHO, PRESS, TEMP = self.hydrostatic_adiabatic_thermo(zfulls)

    THERMODATA = {
      "PRESS": np.tile(PRESS, ntime),
      "TEMP": np.tile(TEMP, ntime),
      "qvap": np.full(int(ngridboxes*ntime), self.qvap),
      "qcond": np.full(int(ngridboxes*ntime), self.qcond), 
      "WVEL": np.array([]), 
      "UVEL": np.array([]),
      "VVEL": np.array([])
    }

    if self.WMAX != None:
      rhotilda = RHO / self.RHO0
      WVEL, UVEL = divfree_flowfield2D(self.WMAX, rhotilda, self.Zlength, 
                                       self.Xlength, zfulls, xfulls)
      THERMODATA["WVEL"] =  np.tile(WVEL, ntime)
      THERMODATA["UVEL"] =  np.tile(UVEL, ntime)

      if self.VVEL != None:
        THERMODATA["VVEL"] =  np.full(int(ngridboxes*ntime), self.VVEL)

    return THERMODATA 