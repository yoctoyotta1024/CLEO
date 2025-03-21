"""
----- CLEO -----
File: pltsds.py
Project: plotssrc
Created Date: Friday 17th November 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
Last Modified: Friday 21st March 2025
Modified By: CB
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
Copyright (c) 2023 MPI-M, Clara Bayley
-----
File Description:
examples for plotting individual superdroplets
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import random
from matplotlib.markers import MarkerStyle

sys.path.append("../../../")  # for imports from pySD package
from pySD.sdmout_src import sdtracing


def maxmin_radii_label(r, allradii):
    """if radius, r, is biggest or smallest
    out of allradii, return appropraite label"""

    label = None
    if r == np.amin(allradii):
        label = "{:.2g}\u03BCm".format(r)
        label = "min r0 = " + label
    elif r == np.amax(allradii):
        label = "{:.2g}\u03BCm".format(r)
        label = "max r0 = " + label

    return label


def individ_radiusgrowths_figure(time, radii, savename=""):
    """plots of droplet radii growth given array of radii
    of shape [time, SDs]"""

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 8))

    for i in range(radii.shape[1]):
        label = maxmin_radii_label(radii[0, i], radii[0, :])
        ax.plot(time, radii[:, i], linewidth=0.8, label=label)

    ax.set_xlabel("time /s")
    ax.set_yscale("log")
    ax.legend(fontsize=13)

    ax.set_ylabel("droplet radius /\u03BCm")

    fig.tight_layout()

    if savename != "":
        fig.savefig(savename, dpi=400, bbox_inches="tight", facecolor="w", format="png")
        print("Figure .png saved as: " + str(savename))

    plt.show()

    return fig, ax


def plot_randomsample_superdrops(time, sddata, totnsupers, nsample, savename=""):
    """plot timeseries of the attributes of a
    random sample of superdroplets"""

    fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(10, 6))
    fig.suptitle("Time Evolution of a Random Sample of Superdroplets")

    minid, maxid = 0, int(totnsupers)  # largest value of ids to sample
    ids2plot = random.sample(list(range(minid, maxid, 1)), nsample)

    attrs = ["radius", "xi", "msol"]
    for a, attr in enumerate(attrs):
        try:
            data = sdtracing.attribute_for_superdroplets_sample(
                sddata, attr, ids=ids2plot
            )
            axs[0, a].plot(time.mins, data, linewidth=0.8)
        except IndexError:
            print("WARNING: didn't plot " + attr)

    mks = MarkerStyle("o", fillstyle="full")
    coords = ["coord3", "coord1", "coord2"]
    for a, coord in enumerate(coords):
        try:
            data = sdtracing.attribute_for_superdroplets_sample(
                sddata, coord, ids=ids2plot
            )
            data = data / 1000  # [km]
            axs[1, a].plot(time.mins, data, linestyle="", marker=mks, markersize=0.2)
        except IndexError:
            print("WARNING: didn't plot " + coord)

    axs[0, 0].set_yscale("log")
    axs[0, 0].set_ylabel("radius, r /\u03BCm")
    axs[0, 1].set_ylabel("multiplicity, \u03BE")
    axs[0, 2].set_ylabel("solute mass, msol /g")
    axs[1, 0].set_ylabel("zcoord /km")
    axs[1, 1].set_ylabel("xcoord /km")
    axs[1, 2].set_ylabel("ycoord /km")
    for ax in axs[1, :]:
        ax.set_xlabel("time /min")

    fig.tight_layout()

    if savename != "":
        fig.savefig(savename, dpi=400, bbox_inches="tight", facecolor="w", format="png")
        print("Figure .png saved as: " + str(savename))

    plt.show()

    return fig, axs


def plot_randomsample_superdrops_2dmotion(
    sddata,
    totnsupers,
    nsample,
    savename="",
    colors=None,
    arrows=False,
    ids2plot=None,
    israndom=True,
    fig=None,
    ax=None,
):
    """plot timeseries of the attributes of a
    random sample of superdroplets"""
    if fig is None:
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6, 6))

    if ids2plot is None:
        minid, maxid = 0, int(totnsupers)  # largest value of ids to sample
        if israndom:
            ids2plot = random.sample(list(range(minid, maxid, 1)), nsample)
        else:
            ids2plot = np.linspace(0, maxid - 1, nsample, dtype=int)

    mks = MarkerStyle("o", fillstyle="full")
    coordz = (
        sdtracing.attribute_for_superdroplets_sample(sddata, "coord3", ids=ids2plot)
        / 1000
    )  # [km]
    coordx = (
        sdtracing.attribute_for_superdroplets_sample(sddata, "coord1", ids=ids2plot)
        / 1000
    )  # [km]

    if colors is not None:
        ax.scatter(coordx, coordz, marker=mks, s=0.4, cmap=colors[0], c=colors[1])
    else:
        ax.plot(coordx, coordz, linestyle="", marker=mks, markersize=0.4)

    if arrows:
        n2plt = min(300, coordx.shape[1])
        drops2arrow = random.sample(list(range(0, coordx.shape[1], 1)), n2plt)
        for n in drops2arrow:  # must loop over drops to get nice positioning of arrows
            x = coordx[:, n][np.logical_not(np.isnan(coordx[:, n]))]
            z = coordz[:, n][np.logical_not(np.isnan(coordz[:, n]))]

            u = np.diff(x)
            w = np.diff(z)
            norm = np.sqrt(u**2 + w**2)
            pos_x = x[:-1] + u / 2
            pos_z = z[:-1] + w / 2

            sl = list(range(0, len(pos_x), 100))
            ax.quiver(
                pos_x[sl],
                pos_z[sl],
                (u / norm)[sl],
                (w / norm)[sl],
                angles="xy",
                zorder=5,
                pivot="mid",
                scale=50,
            )

    ax.set_ylabel("zcoord /km")
    ax.set_xlabel("xcoord /km")

    fig.tight_layout()

    if savename != "":
        fig.savefig(savename, dpi=400, bbox_inches="tight", facecolor="w", format="png")
        print("Figure .png saved as: " + str(savename))

    plt.show()

    return fig, ax


def plot_superdrops_2dmotion(
    sddata,
    ids2plot,
    savename="",
    colors=None,
    arrows=False,
    fig=None,
    ax=None,
):
    fig, ax = plot_randomsample_superdrops_2dmotion(
        sddata,
        np.nan,
        np.nan,
        savename=savename,
        colors=colors,
        arrows=arrows,
        ids2plot=ids2plot,
        israndom=False,
        fig=fig,
        ax=ax,
    )

    return fig, ax
