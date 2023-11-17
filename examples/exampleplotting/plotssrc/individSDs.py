'''
----- CLEO -----
File: individSDs.py
Project: plotssrc
Created Date: Friday 17th November 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
Last Modified: Friday 17th November 2023
Modified By: CB
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
Copyright (c) 2023 MPI-M, Clara Bayley
-----
File Description:
examples for plotting individual superdroplets
'''

import numpy as np
import matplotlib.pyplot as plt

def maxmin_radii_label(r, allradii):
  ''' if radius, r, is biggest or smallest
  out of allradii, return appropraite label'''

  label = None
  if r == np.amin(allradii):
      label = '{:.2g}\u03BCm'.format(r)
      label = 'min r0 = '+label
  elif r == np.amax(allradii):
      label = '{:.2g}\u03BCm'.format(r)
      label = 'max r0 = '+label

  return label

def individ_radiusgrowths_figure(time, radii, savename=""):
  ''' plots of droplet radii growth given array of radii
  of shape [time, SDs] '''
  
  fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 8))
  
  for i in range(radii.shape[1]):
      label = maxmin_radii_label(radii[0,i], radii[0,:])
      ax.plot(time, radii[:,i], linewidth=0.8, label=label)

  ax.set_xlabel('time /s')
  ax.set_yscale('log')
  ax.legend(fontsize=13)

  ax.set_ylabel('droplet radius /\u03BCm')

  fig.tight_layout()

  if savename != "":
    fig.savefig(savename, dpi=400,
                bbox_inches="tight", facecolor='w', format="png")
    print("Figure .png saved as: "+savename)
  
  plt.show()

  return fig, ax


