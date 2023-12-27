'''
----- CLEO -----
File: dryrgens.py
Project: initsuperdropsbinary_src
Created Date: Wednesday 27th December 2023
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
various ways of generating dryradii of
superdroplets for their initial conditions
'''

import numpy as np

class ScaledRadiiGen:
  '''method to generate superdroplet dryradii that
  are the radii divided by a scale factor 'sf' '''
  
  def __init__(self, scale_factor):

      self.sf = scale_factor 
      
  def __call__(self, radii):
      ''' Returns dryradii for nsupers [m]'''

      return radii / self.sf # units [m]