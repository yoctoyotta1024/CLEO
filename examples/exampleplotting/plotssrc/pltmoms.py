"""
----- CLEO -----
File: pltmoms.py
Project: plotssrc
Created Date: Tuesday 21st November 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
Copyright (c) 2023 MPI-M, Clara Bayley
-----
File Description:
examples for ploting momments of the droplet distirbution
"""

import numpy as np
import matplotlib.pyplot as plt


def plot_totnsupers(time, totnsupers, savename="", fig=None, ax=None):
    if fig is None:
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 5), sharex=True)

    ax.plot(time.mins, totnsupers)

    ax.set_ylabel("domain total number of superdroplets")
    ax.set_xlabel("time /min")

    fig.tight_layout()
    if savename != "":
        fig.savefig(savename, dpi=400, bbox_inches="tight", facecolor="w", format="png")
        print("Figure .png saved as: " + str(savename))


def plot_domainmassmoments(time, massmoms, savename=""):
    def totmassmom(massmom):
        """mass moment summed over entire domain"""
        return np.sum(massmom, axis=(1, 2, 3))

    fig, axs = plt.subplots(nrows=5, ncols=1, figsize=(6, 8), sharex=True)
    fig.suptitle("Total Mass Moments Over Domain")

    axs[0].plot(time.mins, totmassmom(massmoms.nsupers))
    axs[1].plot(time.mins, totmassmom(massmoms.mom0))
    axs[2].plot(time.mins, totmassmom(massmoms.mom1))
    axs[3].plot(time.mins, totmassmom(massmoms.mom2))
    meaneffmass = np.mean((massmoms.effmass), axis=(1, 2, 3))
    axs[4].plot(time.mins, meaneffmass)

    axs[0].set_ylabel("number of\nsuperdroplets")
    axs[1].set_ylabel("$\u03BB^{m}_{0}$, number\nof droplets")
    axs[2].set_ylabel("$\u03BB^{m}_{1}$, droplet\nmass /g")
    axs[3].set_ylabel("$\u03BB^{m}_{2}$\n~reflectivity /g$^2$")
    ylab4 = "mean effective\ndroplet mass,\n<$\u03BB^{m}_{2}$/$\u03BB^{m}_{1}>$ /g"
    axs[4].set_ylabel(ylab4)

    axs[-1].set_xlabel("time /min")

    fig.tight_layout()
    if savename != "":
        fig.savefig(savename, dpi=400, bbox_inches="tight", facecolor="w", format="png")
        print("Figure .png saved as: " + str(savename))
