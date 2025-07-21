"""
Copyright (c) 2024 MPI-M, Clara Bayley


----- CLEO -----
File: attrgens_shima2009.py
Project: boxmodelcollisions
Created Date: Saturday 15th June 2024
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
Radius and xi generators as in Shima et al. 2009
"""

import numpy as np


class SampleRadiiShima2009:
    """method to generate superdroplet radii by sampling between rspan[0] and rspan[1]
    with probability weighted by volume exponential probability distribution, as in
    Shima et al. 2009"""

    def __init__(self, radius0, rspan):
        self.rspan = rspan  # [min, max] radii to sample between [m]
        self.vol0 = (
            4.0 / 3.0 * np.pi * radius0**3
        )  # peak of volume exponential distribution [m^3]

    def __call__(self, nsupers):
        """Returns radii for nsupers sampled from rspan [m]"""

        return self.generate_radiisample(nsupers)  # units [m]

    def generate_radiisample(self, nbins):
        """Sample from rspan[0] to rspan[1] with probability of radius[m] given by
        exponential distirubiton in volume"""

        if nbins:
            radii = self.sample_truncated_volume_exponential(nbins)
            print(len(radii))
            return radii  # [m]
        else:
            return np.array([])

    def sample_truncated_volume_exponential(self, nbins):
        """sample volume exponential rejecting samples which result
        in radii outside rspan[0] <= r <= rspan[1] until nbins number of sample
        radii are found"""
        radii_sample = []
        while len(radii_sample) < nbins:
            n = nbins - len(radii_sample)
            vols = np.random.exponential(self.vol0, n)  # [m^3]
            radii = np.power(3.0 / (4.0 * np.pi) * vols, 1.0 / 3.0)
            in_range = np.where(radii >= self.rspan[0], radii, 0)
            in_range = np.where(radii <= self.rspan[1], in_range, 0)
            radii_sample.extend(in_range[in_range > 0])

        return np.array(radii_sample)[:nbins]


class SampleXiShima2009:
    """probability of radius given by
    Shima et al. (2009) is equal for all radii"""

    def __call__(self, radii, totxi):
        return self._probdistrib(radii)

    def _probdistrib(self, radii):
        """Returns probability of each radius in radii according to uniform
        distribution as in Shima et al. 2009"""

        return np.full(radii.shape, 1 / len(radii))
