import numpy as np

from .radiiprobdistribs import *
from ..gbxboundariesbinary_src.read_gbxboundaries import calc_domainvol

class MonoAttrsGen:
    ''' method to generate superdroplets with an
    attribute all equal to attr0 '''

    def __init__(self, attr0):

        self.attr0 = attr0

    def __call__(self, nsupers):
        ''' Returns attribute for nsupers all
        with the value of attr0 '''

        attrs = np.full(nsupers, self.attr0)

        return attrs

class SampleDryradiiGen:
    '''method to generate superdroplet radii by sampling bins bewteen
    rspan that are linearly spaced in log10(r) '''

    def __init__(self, rspan, random):

        self.rspan = rspan
        self.random = random

    def __call__(self, nsupers):
        ''' Returns dryradii for nsupers sampled from rspan [m]'''

        self.nbins = nsupers
        dryradii = self.generate_dryradiisample()

        return dryradii  # units [m]

    def generate_dryradiisample(self, edges=False):
        ''' Divide rspan [m] into evenly spaced bins in log10(r). If edges=True,
        return values of radii at edges of bins. Else sample each bin randomly to
        obtain the dry radius of 'nsupers' no. of superdroplets unless random=False.
        If random=False, return radii evenly distirbuted in log10(r /m) space '''

        if self.nbins:
            log10redgs = np.linspace(np.log10(self.rspan[0]), np.log10(
                                    self.rspan[1]), self.nbins+1)  # log10(r) bin edges

            if edges:
                redgs = 10**(log10redgs)
                return redgs

            else:
                if not self.random:
                    radii = self.centres_log10rbins(log10redgs)

                else:
                    radii = self.randomlysample_log10rbins(log10redgs)

                return radii  # [m]
        else:
            return np.array([])
        
    def centres_log10rbins(self, log10redgs):
        ''' return radii [m] that are at centres of
        evenly spaced bins of log10(radii /m)'''

        log10rcens = 0.5*(log10redgs[1:] + log10redgs[:-1])
        rcens = 10**(log10rcens)

        return rcens  # [m]

    def randomlysample_log10rbins(self, log10redgs):
        ''' given the bin edges, randomly sample evenly spaced bins
        of log10(radii /m) and return the resultant radii [m]'''

        log10r_binwidth = (log10redgs[-1] - log10redgs[0])/self.nbins

        randlog10deltar = np.random.uniform(
            low=0.0, high=log10r_binwidth, size=self.nbins)
        randlog10radii = log10redgs[:-1] + randlog10deltar

        radii = 10**(randlog10radii)

        return radii  # [m]

class MonoCoordGen:
    ''' method to generate superdroplets with 
     coord all equal to coord0 '''

    def __init__(self, coord0):

        self.coord0 = coord0

    def __call__(self, nsupers, coordrange):
        ''' Returns coord for nsupers all
        with the value of coord0 '''

        if (self.coord0 >= coordrange[0] and self.coord0 < coordrange[1]):
            attrs = np.full(nsupers, self.coord0)
        else:
            attrs = np.array([])
        
        return attrs
    
class SampleCoordGen:
    ''' method to generate 'nsupers'
    no. of superdroplets' coord [m]
    by sampling in range bewteen coordspan '''

    def __init__(self, random):

        self.random = random

    def __call__(self, nsupers, coordrange):
        ''' Returns coord3 for nsupers
        sampled from coord3span [m]'''

        if not self.random:
          coord = np.linspace(coordrange[0], coordrange[1],
                              nsupers, endpoint=False)
        else:
          coord = np.random.uniform(low=coordrange[0], high=coordrange[1], 
                                      size=nsupers)

        return coord  # units [m]

class InitManyAttrsGen:
    ''' class for functions to generate attributes of superdroplets
    given the generators for independent attributes
    e.g. for radius and coord3 in substantation of class'''

    def __init__(self, dryradiigen, radiiprobdist,
                 coord3gen, coord1gen, coord2gen):

        self.dryradiigen = dryradiigen
        self.radiiprobdist = radiiprobdist
        
        self.coord3gen = coord3gen
        self.coord1gen = coord1gen
        self.coord2gen = coord2gen

        self.ncoordsgen = sum(x is not None for x in [coord3gen, coord2gen, coord1gen])

    def mass_solutes(self, dryradii, RHO_SOL):
        ''' return the mass [Kg] of the solute in superdroplets given their 
        dry radii [m] and solute density [Kg m^3]'''

        m_sols = 4.0/3.0 * np.pi * (dryradii**3) * RHO_SOL

        return m_sols  # [Kg]

    def multiplicities(self, radii, NUMCONC, samplevol):
        ''' Calculate the multiplicity of the dry radius of each
        superdroplet given it's probability such that the total number
        concentration [m^-3] of real droplets in the volume, vol, [m^3]
        is about 'numconc'. Raise an error if any of the calculated
        mulitiplicities are zero '''

        prob = self.radiiprobdist(radii)
        eps = np.rint(prob * NUMCONC * samplevol)

        if any(eps == 0):
          num = len(eps[eps==0])
          warning = "WARNING, "+str(num)+" out of "+str(len(eps))+" SDs"+\
            " created with multiplicity = 0. Consider increasing numconc"+\
            " or decreaing range of radii sampled."
          raise ValueError(warning)

        return np.array(eps, dtype=np.uint)

    def check_coordsgen_matches_modeldimension(self, SDnspace):
       
        if SDnspace != self.ncoordsgen:
            errmsg = str(self.ncoordsgen)+" coord generators specified "+\
                    "but SDnspace = "+str(SDnspace)
            raise ValueError(errmsg)

    def check_totalnumconc(self, multiplicities, NUMCONC, samplevol):
        '''check number concentration of real droplets calculated from
        multiplicities is same as input value for number conc. Also check
        total number of real droplets is within 10% of the expected
        value given the input number conc and sample volume '''

        nreals = np.rint(NUMCONC * samplevol)
        calcnreals = np.rint(np.sum(multiplicities))
        calcnumconc = np.rint(calcnreals/samplevol)
        
        if np.rint(NUMCONC) != calcnumconc:
          errmsg = "total real droplet concentration"+\
                    " {:0g} != numconc, {:0g}".format(calcnumconc, NUMCONC)
          raise ValueError(errmsg)

        if abs(nreals - calcnreals) > 0.001*nreals:
          errmsg = "total no. real droplets, {:0g},".format(calcnreals)+\
            " not consistent with sample volume {:.3g} m^3".format(samplevol)
          raise ValueError(errmsg)

        else:
          msg = "--- total real droplet concentration = "+\
            "{:0g} m^-3 in {:.3g} m^3 volume --- ".format(calcnumconc, samplevol)
          print(msg)

    def generate_attributes(self, nsupers, RHO_SOL, NUMCONC, gridboxbounds):
        ''' generate superdroplets (SDs) attributes that have dimensions
        by calling the appropraite generating functions'''

        gbxvol = calc_domainvol(gridboxbounds[0:2], gridboxbounds[2:4], 
                                gridboxbounds[4:]) # [m^3]
        
        dryradii = self.dryradiigen(nsupers) # [m]

        mass_solutes = self.mass_solutes(dryradii, RHO_SOL)  # [Kg]

        multiplicities = self.multiplicities(dryradii, NUMCONC, gbxvol)

        if nsupers > 0:  
            self.check_totalnumconc(multiplicities, NUMCONC, gbxvol) 
         
        return multiplicities, dryradii, mass_solutes # units [], [m], [Kg], [m]

    def generate_coords(self, nsupers, SDnspace, gridboxbounds):
        ''' generate superdroplets (SDs) attributes that have dimensions
        by calling the appropraite generating functions'''

        self.check_coordsgen_matches_modeldimension(SDnspace)
       
        coord3, coord1, coord2 = np.array([]), np.array([]), np.array([])
        
        if self.coord3gen:
            coord3range = [gridboxbounds[0], gridboxbounds[1]] # [min,max] coord3 to sample within
            coord3 = self.coord3gen(nsupers, coord3range)

            if self.coord1gen:
                coord1range = [gridboxbounds[2], gridboxbounds[3]] # [min,max] coord1 to sample within
                coord1 = self.coord1gen(nsupers, coord1range)

                if self.coord2gen:
                    coord2range = [gridboxbounds[4], gridboxbounds[5]] # [min,max] coord2 to sample within
                    coord2 = self.coord2gen(nsupers, coord2range)


        return coord3, coord1, coord2 # units [m], [m], [m]