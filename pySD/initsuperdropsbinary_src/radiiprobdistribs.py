import numpy as np

class DiracDelta:
  ''' probability of radius nonzero if it is
  closest value in sample of radii to r0 '''

  def __init__(self, r0):
    
    self.r0 = r0

  def __call__(self, radii):
    ''' Returns probability of radius in radii sample for
    discrete version of dirac delta function centred on 
    value of r in radii closest to r0. For each radius in radii, 
    probability of that radius = 0 if it's not the closest value
    in radii to r0. If it is the closest, the probability is maximal
    (ie. prob = 1 and is then re-normalised such that sum of the
    probalilities over the sample = 1) ''' 

    if radii.any():
      diff = np.abs(radii - self.r0)
      
      probs = np.where(diff == np.min(diff), 1, 0)

      return probs / np.sum(probs) #normalise so sum(prob) = 1
    else:
      return np.array([])
    
class VolExponential:
  ''' probability of radius given by exponential in
  volume distribution as defined by Shima et al. (2009) '''

  def __init__(self, radius0, rspan):
    
    self.radius0 = radius0 # peak of volume exponential distribution [m]
    self.rspan = rspan

  def __call__(self, radii):
    ''' Returns probability of eaach radius in radii according to 
    distribution where probability of volume is exponential and bins 
    for radii are evently spaced in ln(r).
    typical parameter values:
    radius0 = 30.531e-6 # [m]
    numconc = 2**(23) # [m^-3] ''' 

    rwdth = np.log(self.rspan[-1]/self.rspan[0]) # assume equally spaced bins in ln(r)
    
    dn_dV = np.exp(-(radii/self.radius0)**3) # prob density P(volume)

    probs = rwdth * radii**3 * dn_dV # prob density * ln(r) bin width

    return probs / np.sum(probs) # normalise so sum(prob) = 1

class LnNormal:
  ''' probability of radius given by lognormal distribution
  as defined by section 5.2.3 of "An Introduction to clouds from
  the Microscale to Climate" by Lohmann, Luond and Mahrt and radii sampled
  from evenly spaced bins in ln(r).
  typical parameter values:
  geomeans = [0.02e-6, 0.2e-6, 3.5e-6] # [m]               
  geosigs = [1.55, 2.3, 2]                    
  scalefacs = [1e6, 0.3e6, 0.025e6]
  numconc = 1e9 # [m^-3]'''

  def __init__(self, geomeans, geosigs, scalefacs):
    
    nmodes = len(geomeans)
    if nmodes != len(geosigs) or nmodes != len(scalefacs):
      errmsg = "parameters for number of lognormal modes is not consistent"
      raise ValueError(errmsg)
    else:
      self.nmodes = nmodes
      self.geomeans = geomeans
      self.geosigs = geosigs
      self.scalefacs = scalefacs

  def __call__(self, radii):
    ''' Returns probability of each radius in radii derived
    from superposition of Logarithmic (in e) Normal Distributions ''' 

    probs = np.zeros(radii.shape)
    for n in range(self.nmodes):
      probs += self.lnnormaldist(radii, self.scalefacs[n],
                             self.geomeans[n], self.geosigs[n])

    return probs / np.sum(probs) #normalise so sum(prob) = 1 

  def lnnormaldist(self, radii, scalefac, geomean, geosig):
    ''' calculate probability of radii given the paramters of a
    lognormal dsitribution accordin to equation 5.8 of "An 
    Introduction to clouds from the Microscale to Climate" 
    by Lohmann, Luond and Mahrt ''' 

    sigtilda = np.log(geosig)
    mutilda = np.log(geomean)

    norm = scalefac / (np.sqrt(2*np.pi)*sigtilda)
    exponent = -(np.log(radii)-mutilda)**2/(2*sigtilda**2)

    dn_dlnr = norm*np.exp(exponent)                         # eq.5.8 [lohmann intro 2 clouds]
    
    return dn_dlnr
