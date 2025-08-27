"""
----- CLEO -----
File: shima2009fig.py
Project: plotssrc
Created Date: Friday 17th November 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
Copyright (c) 2023 MPI-M, Clara Bayley
-----
File Description:
functions for plotting similar to
figure 2(a) from Shima et al. 2009
"""

import numpy as np
import awkward as ak
import matplotlib.pyplot as plt
from scipy.special import iv

from .pltdist import logr_distribution


def plot_validation_figure(
    witherr,
    time,
    superdrops,
    t2plts,
    domainvol,
    n_a,
    r_a,
    smoothsig,
    xlims=[10, 5000],
    savename="",
    withgol=False,
):
    variables2slice = ["time", "radius", "xi"]
    superdrops_timeslice = superdrops.time_slice(
        t2plts, variables2slice, attach_time=True, time=time.secs, time_units="s"
    )

    nbins = 500
    non_nanradius = ak.nan_to_none(superdrops["radius"])
    rspan = [ak.min(non_nanradius), ak.max(non_nanradius)]

    fig, ax, ax_err = setup_validation_figure(witherr, xlims)

    for n in range(len(superdrops_timeslice.time())):
        t = superdrops_timeslice.time()[n]
        tlab = "t = {:.2f}".format(t) + superdrops_timeslice.time_units()
        c = "C" + str(n)

        if withgol:
            golsol, hcens = golovin_analytical(
                rspan, t, nbins, n_a, r_a, superdrops.RHO_L()
            )
            plot_golovin_analytical_solution(ax, hcens, golsol, n, c)

        radius = superdrops_timeslice["radius"][n]
        xi = superdrops_timeslice["xi"][n]
        hist, hcens = calc_massdens_distrib(
            rspan, nbins, domainvol, xi, radius, superdrops, smoothsig
        )

        if smoothsig:
            ax.plot(hcens, hist, label=tlab, color=c)
        else:
            ax.step(hcens, hist, label=tlab, where="mid", color=c)

        if witherr:
            golsol, hcens = golovin_analytical(
                rspan, t, nbins, n_a, r_a, superdrops.RHO_L()
            )
            diff = hist - golsol
            ax_err.plot(hcens, diff, c=c)

    ax.legend()

    fig.tight_layout()

    if savename != "":
        fig.savefig(savename, dpi=400, bbox_inches="tight", facecolor="w", format="png")
        print("Figure .png saved as: " + str(savename))

    return fig, ax


def setup_validation_figure(witherr, xlims):
    if witherr:
        gd = dict(height_ratios=[5, 1])
        fig, [ax, ax_err] = plt.subplots(
            ncols=1, nrows=2, figsize=(8, 7), gridspec_kw=gd
        )
    else:
        fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(8, 6))

    ax.set_xscale("log")
    ax.set_xlim(xlims)
    ax.set_xlabel("radius, r, /\u03BCm")

    ax.set_ylabel("mass density distribution,\n g(lnR) /g m$^{-3}$ / unit lnR")
    ax.set_ylim([0, 1.8])

    if witherr:
        ax_err.set_xscale("log")
        ax_err.set_xlim(xlims)
        ax_err.set_xlabel("radius, r, /\u03BCm")
        ax.set_xlabel(None)
        ax.set_xticklabels([])

        ax_err.set_ylabel("error /g m$^{-3}$\n/ unit lnR")
        ax_err.set_ylim([-0.25, 0.25])

        return fig, ax, ax_err

    else:
        return fig, ax, []


def golovin_analytical(rspan, t, nbins, n_a, r_a, RHO_L):
    b = 1500
    rspan = [r / 1e6 for r in rspan]  # convert from microns to m

    if t < 1e-16:
        t = 1e-10

    hedgs = np.linspace(
        np.log(rspan[0]), np.log(rspan[1]), nbins + 1
    )  # edges to lnr bins
    hcens = (hedgs[1:] + hedgs[:-1]) / 2  # lnr bin centres
    r = np.exp(hcens)

    vol_a = 4 / 3 * np.pi * r_a**3
    x = 4 / 3 * np.pi * r**3 / vol_a
    tau = 1 - np.exp(-b * n_a * vol_a * t)
    bsl_exp = iv(1, 2 * x * np.sqrt(tau)) * np.exp(-(1 + tau) * x)
    asym = 1 / (2 * np.sqrt(np.pi * x)) * np.exp(x * (2 * np.sqrt(tau) - 1 - tau))
    np.nan_to_num(bsl_exp, copy=False, nan=asym, posinf=asym, neginf=np.inf)

    phi = (1 - tau) / (x * np.sqrt(tau)) * bsl_exp
    n = n_a / vol_a * phi
    # dv_dlnR = 4/3*np.pi*((rwdths[1:])**3-(rwdths[:-1])**3)/hwdths
    dv_dlnR = 3 * (4 / 3 * np.pi * r**3)
    massdens = (
        n * RHO_L * (4 / 3 * np.pi * r**3) * dv_dlnR
    )  # mass density as if water [Kg m^-3 /unit lnR]]

    return massdens * 1000, r * 1e6  # units: [g m^-3 /unit lnR] , [microns]


def plot_golovin_analytical_solution(ax, hcens, golsol, n, c):
    if n == 0:
        # add legend to analytical solution
        glab = "analytical solution"
        ax.plot(hcens, golsol, label=glab, color="k", linestyle="--")

    ax.plot(hcens, golsol, label=None, color=c, linestyle="--")

    return ax


def calc_massdens_distrib(rspan, nbins, domainvol, xi, radius, superdrops, smoothsig):
    m_asif_water = (
        superdrops.vol(radius) * superdrops.RHO_L()
    )  # superdrops mass as if water [g]
    weights = xi * m_asif_water * 1000 / domainvol  # real droplets [g/m^3]

    hist, hedges, hcens = logr_distribution(
        rspan, nbins, radius, weights, perlogR=True, smooth=smoothsig
    )

    return hist, hcens
