'''
----- CLEO -----
File: pltmoms.py
Project: plotssrc
Created Date: Tuesday 21st November 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
Last Modified: Tuesday 21st November 2023
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



