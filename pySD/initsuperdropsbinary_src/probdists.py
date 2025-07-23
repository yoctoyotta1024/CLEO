"""
----- CLEO -----
File: probdists.py
Project: initsuperdropsbinary_src
Created Date: Wednesday 22nd November 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
Copyright (c) 2023 MPI-M, Clara Bayley
-----
File Description:
Class calls' return normalised probability
of radii for various probability distributions
assuming bins are evenly spaced in log10(r)
"""

import numpy as np
from scipy import special


class MinXiDistrib:
    """probability of radius for a given probability distribution but with a
    minimum value such that xi>='xi_min'"""

    def __init__(self, probdistrib, xi_min):
        self.probdistrib = probdistrib
        self.xi_min = xi_min

    def __call__(self, radii, totxi):
        """returns probability for radii from a certain distribution or the
        probability such that xi has minimum value 'xi_min'"""
        prob = self.probdistrib(radii, totxi)
        prob_min = self.xi_min / totxi
        return np.where(prob < prob_min, prob_min, prob)


class CombinedRadiiProbDistribs:
    """probability of radius from the sum of several
    probability distributions"""

    def __init__(self, probdistribs, scalefacs):
        self.probdistribs = probdistribs
        self.scalefacs = scalefacs

        if len(scalefacs) != len(probdistribs):
            errmsg = "relative height of each probability distribution must be given"
            raise ValueError(errmsg)

    def __call__(self, radii, totxi):
        return self._probdistrib(radii, totxi)

    def _probdistrib(self, radii, totxi):
        """returns distribution for radii given by the
        sum of the distributions in probdistribs list"""

        probs = np.zeros(radii.shape)
        for distrib, sf in zip(self.probdistribs, self.scalefacs):
            probs += sf * distrib(radii, totxi)

        return probs / np.sum(probs)  # normalise so sum(prob) = 1


class DiracDelta:
    """probability of radius nonzero if it is
    closest value in sample of radii to r0"""

    def __init__(self, r0):
        self.r0 = r0

    def __call__(self, radii, totxi):
        return self._probdistrib(radii)

    def _probdistrib(self, radii):
        """Returns probability of radius in radii sample for
        discrete version of dirac delta function centred on
        value of r in radii closest to r0. For each radius in radii,
        probability of that radius = 0 if it's not the closest value
        in radii to r0. If it is the closest, the probability is maximal
        (ie. prob = 1 and is then re-normalised such that sum of the
        probalilities over the sample = 1)"""

        if radii.any():
            diff = np.abs(radii - self.r0)

            probs = np.where(diff == np.min(diff), 1, 0)

            return probs / np.sum(probs)  # normalise so sum(prob) = 1
        else:
            return np.array([])


class VolExponential:
    """probability of radius given by exponential in
    volume distribution as defined by Shima et al. (2009)"""

    def __init__(self, radius0, rspan):
        self.radius0 = radius0  # peak of volume exponential distribution [m]
        self.rspan = rspan

    def __call__(self, radii, totxi):
        return self._probdistrib(radii)

    def _probdistrib(self, radii):
        """Returns probability of eaach radius in radii according to
        distribution where probability of volume is exponential and bins
        for radii are evently spaced in ln(r).
        typical parameter values:
        radius0 = 30.531e-6 # [m]
        numconc = 2**(23) # [m^-3]"""

        rwdth = np.log(
            self.rspan[-1] / self.rspan[0]
        )  # assume equally spaced bins in ln(r)

        dn_dV = np.exp(-((radii / self.radius0) ** 3))  # prob density P(volume)

        probs = rwdth * radii**3 * dn_dV  # prob density * ln(r) bin width

        return probs / np.sum(probs)  # normalise so sum(prob) = 1


class LnNormal:
    """probability of radius given by lognormal distribution
    as defined by section 5.2.3 of "An Introduction to clouds from
    the Microscale to Climate" by Lohmann, Luond and Mahrt and radii
    sampled from evenly spaced bins in ln(r).
    typical parameter values:
    geomeans = [0.02e-6, 0.2e-6, 3.5e-6] # [m]
    geosigs = [1.55, 2.3, 2]
    scalefacs = [1e6, 0.3e6, 0.025e6]
    numconc = 1e9 # [m^-3]"""

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

    def __call__(self, radii, totxi):
        return self._probdistrib(radii)

    def _probdistrib(self, radii):
        """Returns probability of each radius in radii derived
        from superposition of Logarithmic (in e) Normal Distributions"""

        probs = np.zeros(radii.shape)
        for n in range(self.nmodes):
            probs += self.lnnormaldist(
                radii, self.scalefacs[n], self.geomeans[n], self.geosigs[n]
            )

        return probs / np.sum(probs)  # normalise so sum(prob) = 1

    def lnnormaldist(self, radii, scalefac, geomean, geosig):
        """calculate probability of radii given the paramters of a
        lognormal dsitribution accordin to equation 5.8 of "An
        Introduction to clouds from the Microscale to Climate"
        by Lohmann, Luond and Mahrt"""

        sigtilda = np.log(geosig)
        mutilda = np.log(geomean)

        norm = scalefac / (np.sqrt(2 * np.pi) * sigtilda)
        exponent = -((np.log(radii) - mutilda) ** 2) / (2 * sigtilda**2)

        dn_dlnr = norm * np.exp(exponent)  # eq.5.8 [lohmann intro 2 clouds]

        return dn_dlnr


class ClouddropsHansenGamma:
    """probability of radius according to gamma distribution for
    shallow cumuli cloud droplets from Poertge et al. 2023"""

    def __init__(self, reff, nueff):
        self.reff = reff
        self.nueff = nueff

    def __call__(self, radii, totxi):
        return self._probdistrib(radii)

    def _probdistrib(self, radii):
        """return gamma distribution for cloud droplets
        given radius [m] using parameters from Poertge
        et al. 2023 for shallow cumuli (figure 12).
        typical values:
        reff = 7e-6 #[m]
        nueff = 0.08 # []"""

        xp = (1 - 2 * self.nueff) / self.nueff
        n0const = (self.reff * self.nueff) ** (-xp)
        n0const = n0const / special.gamma(xp)

        term1 = radii ** ((1 - 3 * self.nueff) / self.nueff)
        term2 = np.exp(-radii / (self.reff * self.nueff))

        probs = n0const * term1 * term2  # dn_dr [prob m^-1]

        return probs / np.sum(probs)  # normalise so sum(prob) = 1


class RaindropsGeoffroyGamma:
    """probability of radius given gamma distribution for
    shallow cumuli rain droplets from Geoffroy et al. 2014"""

    def __init__(self, nrain, qrain, dvol):
        self.nrain = nrain  # raindrop concentration [ndrops/m^3]
        self.qrain = qrain  # rainwater content [g/m^3]
        self.dvol = dvol  # volume mean raindrop diameter [m]

    def __call__(self, radii, totxi):
        return self._probdistrib(radii)

    def _probdistrib(self, radii):
        """returns probability of each radius according to a
        gamma distribution for rain droplets using parameters
        from Geoffroy et al. 2014 for precipitating shallow
        cumuli RICO (see figure 3 and equations 2,3 and 5).
        typical parameter values:
        nrain = 3 / 0.001 # [ndrops/m^3]
        qrain = 0.9 # [g/m^3]
        dvol = 800e-6 #[m]"""

        nu = 18 / ((self.nrain * self.qrain) ** 0.25)  # []
        lamda = (nu * (nu + 1) * (nu + 2)) ** (1 / 3) / self.dvol  # [m^-1]
        const = self.nrain * lamda**nu / special.gamma(nu)

        diam = 2 * radii  # [m]
        probs = const * diam ** (nu - 1) * np.exp(-lamda * diam)  # dn_dr [prob m^-1]

        return probs / np.sum(probs)  # normalise so sum(prob) = 1
