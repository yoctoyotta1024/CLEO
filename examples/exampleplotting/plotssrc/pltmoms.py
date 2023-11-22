'''
----- CLEO -----
File: pltmoms.py
Project: plotssrc
Created Date: Tuesday 21st November 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
Last Modified: Wednesday 22nd November 2023
Modified By: CB
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
Copyright (c) 2023 MPI-M, Clara Bayley
-----
File Description:
examples for ploting momments of the droplet distirbution
'''

import numpy as np
import matplotlib.pyplot as plt

def plot_totnsupers(time, totnsupers, savename=""):
  
  fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8,5), sharex=True)
    
  ax.plot(time.mins, totnsupers)
  
  ax.set_ylabel("domain total number of superdroplets")
  ax.set_xlabel("time /min")
    
  fig.tight_layout()
  if savename != "":
    fig.savefig(savename, dpi=400, bbox_inches="tight",
                facecolor='w', format="png")
    print("Figure .png saved as: "+savename)
  
  plt.show()

def plot_domainmassmoments(time, massmoms, nsupers=False, savename=""):
  ''' TODO '''

  if nsupers != False:
    fig, axs = plt.subplots(nrows=5, ncols=1, figsize=(10,6), sharex=True)
  else:
    fig, axs = plt.subplots(nrows=4, ncols=1, figsize=(8,6), sharex=True)



  
  def totmassmom(massmom):
    '''mass moment summed over entire domain'''
    return  np.sum(massmom, axis=(1,2,3))
  
  fig, axs = plt.subplots(nrows=2, ncols=2, figsize=figsize,
                          sharex=True)
  fig.suptitle("Mass Moments Over Entire Domain")
  
  if nsupers != False:
    print("plot nsueprs too")

  l0 = axs[0].plot(time, totmassmom(massmoms.mom0))
  l1 = axs[1].plot(time, totmassmom(massmoms.mom1)) 
  l2 = axs[2].plot(time, totmassmom(massmoms.mom2))
  meaneffmass = np.mean((massmoms.effmass), axis=(1,2,3))
  l3 = axs[3].plot(time, meaneffmass)

  axs[0].set_ylabel("$\u03BB^{m}_{0}$, number of  droplets")
  axs[1].set_ylabel("$\u03BB^{m}_{1}$, droplet mass /g")
  axs[2].set_ylabel("$\u03BB^{m}_{2}$ ~reflectivity /g$^2$")
  ylab3 = "mean effective droplet mass,\n<$\u03BB^{m}_{2}$/$\u03BB^{m}_{1}>$ /g"
  axs[3].set_ylabel(ylab3)
  
  axs[3].set_xlabel("time /min")

  fig.tight_layout()
  if savename != "":
    fig.savefig(savename, dpi=400, bbox_inches="tight",
                facecolor='w', format="png")
    print("Figure .png saved as: "+savename)
  
  plt.show()
