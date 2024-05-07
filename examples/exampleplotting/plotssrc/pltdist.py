"""
----- CLEO -----
File: pltdist.py
Project: plotssrc
Created Date: Thursday 23rd October 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
Last Modified: Sunday 26th November 2023
Modified By: CB
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
Copyright (c) 2023 MPI-M, Clara Bayley
-----
File Description:
mass, number concentration etc. plots against
radius in evenly spaced ln(radius) bins
"""

import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append("../../../")  # for imports from pySD package
from pySD.sdmout_src import sdtracing


def gaussian_kernel_smoothing(hist, hcens, sig):
    smoothhist = []
    for h in range(len(hist)):
        kernel = (
            1
            / (np.sqrt(2 * np.pi) * sig)
            * np.exp(-((hcens - hcens[h]) ** 2) / (2 * sig**2))
        )
        kernel = kernel / np.sum(kernel)
        smoothhist.append(np.sum(hist * kernel))

    smoothhist = np.asarray(smoothhist)
    smoothhist = np.where(smoothhist < 1e-16, 0, smoothhist)

    return smoothhist, hcens


def logr_distribution(rspan, nbins, radius, wghts, perlogR=False, smooth=False):
    """get distribution of data with weights 'wghts' against
    logr. Uses np.histogram to get frequency of a particular
    value of data that falls in each ln(r) -> ln(r) + dln(r) bin.
    Apply gaussian kernel smoothing if wanted. Note log base e not 10!"""

    # create lnr bins (linearly spaced in lnr)
    hedgs = np.linspace(
        np.log(rspan[0]), np.log(rspan[1]), nbins + 1
    )  # edges to lnr bins
    logrwdths = hedgs[1:] - hedgs[:-1]  # lnr bin widths
    hcens = np.log((np.exp(hedgs[1:]) + np.exp(hedgs[:-1])) / 2)  # lnr bin centres

    # get number frequency in each bin
    hist, hedgs = np.histogram(np.log(radius), bins=hedgs, weights=wghts, density=None)

    if perlogR is True:  # get frequency / bin width
        hist = hist / logrwdths

    if smooth:
        hist, hcens = gaussian_kernel_smoothing(hist, hcens, smooth)

    return hist, np.exp(hedgs), np.exp(hcens)  # units of hedgs and hcens [microns]


def massdens_distrib(radius, xi, mass, vol, rspan, nbins, perlogR, smooth):
    weights = xi * mass / vol  # real droplets [g/m^3]

    hist, hedges, hcens = logr_distribution(
        rspan, nbins, radius, weights, perlogR=perlogR, smooth=smooth
    )

    return hcens, hist


def nsupers_distrib(radius, xi, mass, vol, rspan, nbins, perlogR, smooth):
    weights = None  # number superdroplets []
    hist, hedges, hcens = logr_distribution(
        rspan, nbins, radius, weights, perlogR=perlogR, smooth=smooth
    )

    return hcens, hist


def numconc_distrib(radius, xi, mass, vol, rspan, nbins, perlogR, smooth):
    weights = xi / vol / 1e6  # real droplets [/cm^3]

    hist, hedges, hcens = logr_distribution(
        rspan, nbins, radius, weights, perlogR=perlogR, smooth=smooth
    )

    return hcens, hist


def plot_dists(
    ax,
    distribcalc,
    timesecs,
    data2plt,
    t2plts,
    vol,
    rspan,
    nbins,
    masscalc=None,
    smoothsig=False,
    perlogR=True,
    ylog=False,
):
    for t, tplt in enumerate(t2plts):
        ind = np.argmin(abs(timesecs - tplt))
        tlab = "t = {:.2f}s".format(timesecs[ind])
        c = "C" + str(t)

        radius = data2plt["radius"][t]
        xi = data2plt["xi"][t]
        if masscalc:
            msol = data2plt["msol"][t]
            mass = masscalc(radius, msol)
        else:
            mass = None

        hcens, hist = distribcalc(
            radius, xi, mass, vol, rspan, nbins, perlogR=perlogR, smooth=smoothsig
        )
        if smoothsig:
            ax.plot(hcens, hist, label=tlab, color=c)
        else:
            ax.step(hcens, hist, label=tlab, where="mid", color=c)

    if ylog:
        ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_xlabel("radius, r, /\u03BCm")
    ax.legend()

    return ax


def plot_domainmass_distribs(
    timesecs,
    sddata,
    t2plts,
    domainvol,
    rspan,
    nbins,
    smoothsig=False,
    perlogR=True,
    ylog=False,
    savename="",
):
    fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(8, 6))

    attrs2sel = ["radius", "xi", "msol"]
    data2plt = sdtracing.attributes_at_times(sddata, timesecs, t2plts, attrs2sel)

    plot_dists(
        ax,
        massdens_distrib,
        timesecs,
        data2plt,
        t2plts,
        domainvol,
        rspan,
        nbins,
        sddata.mass,
        smoothsig=smoothsig,
        perlogR=perlogR,
        ylog=ylog,
    )

    if perlogR:
        ax.set_ylabel("droplet mass distribution,\n g(lnR) /g m$^{-3}$ / unit lnR")
    else:
        ax.set_ylabel("droplet mass distribution,\n M(lnR) /g m$^{-3}$")

    fig.tight_layout()

    if savename != "":
        fig.savefig(savename, dpi=400, bbox_inches="tight", facecolor="w", format="png")
        print("Figure .png saved as: " + savename)

    plt.show()

    return fig, ax


def plot_domainnsupers_distribs(
    timesecs,
    sddata,
    t2plts,
    domainvol,
    rspan,
    nbins,
    smoothsig=False,
    perlogR=True,
    ylog=False,
    savename="",
):
    fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(8, 6))

    attrs2sel = ["radius", "xi", "msol"]
    data2plt = sdtracing.attributes_at_times(sddata, timesecs, t2plts, attrs2sel)

    plot_dists(
        ax,
        nsupers_distrib,
        timesecs,
        data2plt,
        t2plts,
        domainvol,
        rspan,
        nbins,
        None,
        smoothsig=smoothsig,
        perlogR=perlogR,
        ylog=ylog,
    )

    if perlogR:
        ax.set_ylabel("number of superdroplets / unit lnR")
    else:
        ax.set_ylabel("number of superdroplets")

    fig.tight_layout()

    if savename != "":
        fig.savefig(savename, dpi=400, bbox_inches="tight", facecolor="w", format="png")
        print("Figure .png saved as: " + savename)

    plt.show()

    return fig, ax


def plot_domainnumconc_distribs(
    timesecs,
    sddata,
    t2plts,
    domainvol,
    rspan,
    nbins,
    smoothsig=False,
    perlogR=True,
    ylog=False,
    savename="",
):
    fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(8, 6))

    attrs2sel = ["radius", "xi", "msol"]
    data2plt = sdtracing.attributes_at_times(sddata, timesecs, t2plts, attrs2sel)

    plot_dists(
        ax,
        numconc_distrib,
        timesecs,
        data2plt,
        t2plts,
        domainvol,
        rspan,
        nbins,
        None,
        smoothsig=smoothsig,
        perlogR=perlogR,
        ylog=ylog,
    )

    if perlogR:
        ax.set_ylabel("real droplet concentration /cm$^{-3}$ / unit lnR")
    else:
        ax.set_ylabel("real droplet concentration /cm$^{-3}$")

    fig.tight_layout()

    if savename != "":
        fig.savefig(savename, dpi=400, bbox_inches="tight", facecolor="w", format="png")
        print("Figure .png saved as: " + savename)

    plt.show()

    return fig, ax
