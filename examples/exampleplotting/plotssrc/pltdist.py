
import numpy as np
import matplotlib.pyplot as plt

def gaussian_kernel_smoothing(hist, hcens, sig):

    smoothhist = []
    for h in range(len(hist)):
        kernel = 1/(np.sqrt(2*np.pi)*sig)*np.exp(-(hcens - hcens[h])**2/(2*sig**2))
        kernel = kernel/np.sum(kernel)
        smoothhist.append(np.sum(hist*kernel))  

    smoothhist = np.asarray(smoothhist)
    smoothhist = np.where(smoothhist<1e-16, 0, smoothhist)  

    return smoothhist, hcens

def logr_distribution(rspan, nbins, radius, wghts,
                      perlogR=False, smooth=False):
  ''' get distribution of data with weights 'wghts' against 
  logr. Uses np.histogram to get frequency of a particular
  value of data that falls in each ln(r) -> ln(r) + dln(r) bin.
  Apply gaussian kernel smoothing if wanted. Note log base e not 10! '''

  # create lnr bins (linearly spaced in lnr)
  hedgs = np.linspace(np.log(rspan[0]), np.log(rspan[1]), nbins+1)  # edges to lnr bins
  logrwdths = hedgs[1:]- hedgs[:-1]                                 # lnr bin widths
  hcens = np.log((np.exp(hedgs[1:])+np.exp(hedgs[:-1]))/2)          # lnr bin centres

  # get number frequency in each bin
  hist, hedgs = np.histogram(np.log(radius), bins=hedgs, 
                              weights=wghts, density=None)
  
  if perlogR == True: # get frequency / bin width
      hist = hist/logrwdths

  if smooth:
    hist, hcens = gaussian_kernel_smoothing(hist, hcens, smooth)

  return hist, np.exp(hedgs), np.exp(hcens) # units of hedgs and hcens [microns]