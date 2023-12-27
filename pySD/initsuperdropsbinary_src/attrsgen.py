'''
----- CLEO -----
File: attrsgen.py
Project: initsuperdropsbinary_src
Created Date: Friday 13th October 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
Last Modified: Wednesday 27th December 2023
Modified By: CB
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
Copyright (c) 2023 MPI-M, Clara Bayley
-----
File Description:
attrsgen generates multiple superdroplet
attributes given individual generators
'''

import numpy as np

from ..gbxboundariesbinary_src import read_gbxboundaries as rgrid

class AttrsGenerator:
    ''' class for functions to generate attributes of
    superdroplets given the generators for independent
    attributes e.g. for radius and coord3 in substantation
    of class'''

    def __init__(self, radiigen, dryradiigen, xiprobdist, 
                 coord3gen, coord1gen, coord2gen):

        self.radiigen = radiigen            # generates radius (solute + water)
        self.dryradiigen = dryradiigen      # generates dry radius (-> solute mass) 
        self.xiprobdist = xiprobdist        # droplet size distribution (-> multiplicity)
        
        self.coord3gen = coord3gen          # generates spatial coordinate 
        self.coord1gen = coord1gen
        self.coord2gen = coord2gen

        self.ncoordsgen = sum(x is not None for x in
                               [coord3gen, coord2gen, coord1gen])

    def mass_solutes(self, nsupers, radii, RHO_SOL):
        ''' return the mass [Kg] of the solute in superdroplets given their 
        dry radii [m] and solute density [Kg m^3]'''

        dryradii = self.dryradiigen(nsupers) # [m]
        dryradii = np.where(radii < dryradii, radii, dryradii)
        msols = 4.0/3.0 * np.pi * (dryradii**3) * RHO_SOL

        return msols  # [Kg]

    def multiplicities(self, radii, NUMCONC, samplevol):
        ''' Calculate the multiplicity of the dry radius of each
        superdroplet given it's probability such that the total number
        concentration [m^-3] of real droplets in the volume, vol, [m^3]
        is about 'numconc'. Raise an error if any of the calculated
        mulitiplicities are zero '''

        prob = self.xiprobdist(radii) # normalised prob distrib
        xi = np.rint(prob * NUMCONC * samplevol)

        if any(xi == 0):
          num = len(xi[xi==0])
          errmsg = "ERROR, "+str(num)+" out of "+str(len(xi))+" SDs"+\
            " created with multiplicity = 0. Consider increasing numconc"+\
            " or changing range of radii sampled."
          raise ValueError(errmsg)

        return np.array(xi, dtype=np.uint)

    def check_coordsgen_matches_modeldimension(self, nspacedims):
       
        if nspacedims != self.ncoordsgen:
            errmsg = str(self.ncoordsgen)+" coord generators specified "+\
                    "but nspacedims = "+str(nspacedims)
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

        gbxvol = rgrid.calc_domainvol(gridboxbounds[0:2],
                                      gridboxbounds[2:4], 
                                      gridboxbounds[4:]) # [m^3]
        
        radii = self.radiigen(nsupers) # [m]

        mass_solutes = self.mass_solutes(nsupers, radii, RHO_SOL)  # [Kg]

        multiplicities = self.multiplicities(radii, NUMCONC, gbxvol)

        if nsupers > 0:  
            self.check_totalnumconc(multiplicities, NUMCONC, gbxvol) 
         
        return multiplicities, radii, mass_solutes # units [], [m], [Kg], [m]

    def generate_coords(self, nsupers, nspacedims, gridboxbounds):
        ''' generate superdroplets (SDs) attributes that have dimensions
        by calling the appropraite generating functions'''

        self.check_coordsgen_matches_modeldimension(nspacedims)
       
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