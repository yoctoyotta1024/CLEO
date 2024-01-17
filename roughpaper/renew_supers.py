'''
----- CLEO -----
File: renew_supers.py
Project: roughpaper
Created Date: Monday 18th December 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
Last Modified: Thursday 28th December 2023
Modified By: CB
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
Copyright (c) 2023 MPI-M, Clara Bayley
-----
File Description:
'''

import numpy as np
from scipy import special
import matplotlib.pyplot as plt

def Log10Redges(rspan, nbins):

  return np.linspace(np.log10(rspan[0]),
                                   np.log10(rspan[1]),
                                   nbins+1)  # log10(r) bin edges

def random_gbxindex(gbxindexes):
  ''' returns [lower, upper] limits of
  randomly selected bin uniformly
  distributed in log10(r) space '''

  i = np.random.randint(0, len(gbxindexes))

  return gbxindexes[i]

def random_radiusbin(log10redges):
  ''' returns [lower, upper] limits of
  randomly selected bin uniformly
  distributed in log10(r) space '''

  i = np.random.randint(0, len(log10redges)-1)

  return log10redges[i], log10redges[i+1]

def sample_radius(log10rlow, log10rup):

  frac = np.random.rand(1)
  log10r = log10rlow + frac * (log10rup - log10rlow)

  return 10**log10r

def sample_coord3(gbxindex, coord3bounds):

  llim, ulim = coord3bounds[gbxindex]
  frac = np.random.rand(1)
  coord3 = llim + frac * (ulim-llim)

  return coord3

class SampleXi:

  def __init__(self, prob_distrib):

    self.prob_distrib = prob_distrib

  def __call__(self, log10rlow, log10rup, numconc, vol):

    rlow, rup = 10**log10rlow, 10**log10rup
    deltar = rup - rlow

    log10r = 0.5*(log10rup + log10rlow)
    radius = 10**log10r

    prob = self.prob_distrib(radius) * deltar

    xi = prob * numconc * vol

    return xi

class ProbDistrib:

  def __init__(self):

    self.reff = 7e-6                     # effective radius [m]
    self.nueff = 0.08                    # effective variance

  def normalised_prob_distrib(self, radius):
    ''' return gamma distribution for cloud droplets given
    radius [m] using parameters from Poertge et al. 2023
    for shallow cumuli (figure 12). normalised such that
    integral of probs over dr = 1.
    typical values:
    reff = 7e-6 #[m]
    nueff = 0.08 # [] '''

    xp = (1-2*self.nueff)/self.nueff
    n0const = (self.reff*self.nueff)**(-xp)
    n0const = n0const / special.gamma(xp)

    term1 = radius**((1-3*self.nueff)/self.nueff)
    term2 = np.exp(-radius/(self.reff*self.nueff))

    probs = n0const * term1 * term2 # dn_dr [prob m^-1]

    return probs # normalised probability of droplet in range r -> r+dr

  def __call__(self, radius):

    return self.normalised_prob_distrib(radius)

class RenewSuperdrop:

  def __init__(self, gbxindexes, rspan, nbins, numconc):

    self.gbxindexes = gbxindexes
    self.log10redges = Log10Redges(rspan, nbins)
    self.numconc = numconc

  def __call__(self, gbxmaps):

    coord3bounds = gbxmaps[0]
    gbxvols = gbxmaps[1]

    gbxindex = random_gbxindex(self.gbxindexes)
    coord3 = sample_coord3(gbxindex, coord3bounds)

    log10rlow, log10rup = random_radiusbin(self.log10redges)
    radius = sample_radius(log10rlow, log10rup)
    xi = sample_xi(log10rlow, log10rup, numconc, gbxvols[gbxindex])

    return Superdroplet(xi, radius, gbxindex, coord3)

##### ------------------------------------------------------------ #####

class MoveSupers:

  def __init__(self, gbxmaps, renew_superdrop):

    self.renew_superdrop = renew_superdrop
    self.gbxmaps = gbxmaps

  def backwards_coord3idx(self, at_cartesiandomainboundary, drop):

    if (at_cartesiandomainboundary):

      drop = self.renew_superdrop(self.gbxmaps)

    return drop


class Superdroplet:

  def __init__(self, xi, radius, gbxindex, coord3):

    self.xi = xi
    self.radius = radius
    self.gbxindex = gbxindex
    self.coord3 = coord3

##### ------------------------------------------------------------ #####

nsupers = 100
ngbxs = 5
rspan = [2e-7, 3e-5] # [m]
numconc = 100 * 1e6 # 100/cm^3

gbxindexes = range(0, ngbxs)
log10redges = Log10Redges(rspan, nsupers)
sample_xi = SampleXi(ProbDistrib())

coord3bounds = {}
gbxvols = {}
lim1 = 0.0
lim2 = 100.0
for gbxindex in gbxindexes:
  coord3bounds[gbxindex] = [lim1, lim2]
  gbxvols[gbxindex] = (lim2 - lim1) * 50 * 50
  lim1 = lim2
  lim2 += 100
  print("gbx: ", gbxindex, coord3bounds[gbxindex], gbxvols[gbxindex])

renew_superdrop = RenewSuperdrop(gbxindexes, rspan, nsupers, numconc)
gbxmaps = [coord3bounds, gbxvols]
movesupers = MoveSupers(gbxmaps, renew_superdrop)
superdrops = []
for i in range(nsupers):
  drop = Superdroplet(0, 0, 0, 0)
  superdrops.append(movesupers.backwards_coord3idx(True, drop))

##### ------------------------------------------------------------ #####

fig, axs = plt.subplots(nrows=2, ncols=4)
axs = axs.flatten()
sd_gbxindex = []
sd_coord3 = []
sd_radius = []
sd_xi = []
for i in range(nsupers):
  sd_gbxindex.append(superdrops[i].gbxindex)
  sd_coord3.append(superdrops[i].coord3)
  sd_radius.append(superdrops[i].radius)
  sd_xi.append(superdrops[i].xi)
sd_gbxindex = np.asarray(sd_gbxindex).flatten()
sd_coord3 = np.asarray(sd_coord3).flatten()
sd_radius = np.asarray(sd_radius).flatten()
sd_xi = np.asarray(sd_xi).flatten()

axs[0].scatter(sd_gbxindex, sd_gbxindex)
axs[1].scatter(sd_gbxindex, sd_coord3)
axs[2].scatter(sd_gbxindex, sd_radius)
axs[3].scatter(sd_radius, sd_xi)

hedges = np.arange(0, ngbxs+1, 1)
hist, hedges = np.histogram(sd_gbxindex, bins=hedges)
hcens = (hedges[1:]+hedges[:-1])/2
print("mean nsupers per gbx:", np.mean(hist))
axs[4].plot(hcens, hist)

hedges = np.arange(0, lim2, 100)
hist, hedges = np.histogram(sd_coord3, bins=hedges)
hcens = (hedges[1:]+hedges[:-1])/2
print(hedges, hcens)
print("mean nsupers per bin:", np.mean(hist))
axs[5].plot(hcens, hist)

hedges = 10**Log10Redges(rspan, 25)
hist, hedges = np.histogram(sd_radius, bins=hedges)
hcens = (hedges[1:]+hedges[:-1])/2
axs[6].plot(hcens, hist)

hedges = 10**Log10Redges(rspan, 25)
hist, hedges = np.histogram(sd_radius, bins=hedges, weights=sd_xi)
hcens = (hedges[1:]+hedges[:-1])/2
axs[7].plot(hcens, hist)

print("sum: ", np.sum(sd_xi), np.sum(sd_xi)/gbxvols[0]/1e6)
axs[2].set_yscale("log")
axs[3].set_xscale("log")
axs[6].set_xscale("log")
axs[7].set_xscale("log")
plt.savefig("./sample.png")
plt.show()
