"""
----- CLEO -----
File: pltsds.py
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
examples for plotting individual superdroplets
"""

import awkward as ak
import numpy as np
import matplotlib.pyplot as plt
import random
import warnings
from matplotlib.markers import MarkerStyle


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

    radii0 = [r[0] for r in radii]
    for i in range(len(radii)):
        label = maxmin_radii_label(radii[i][0], radii0)
        ax.plot(time, radii[i], linewidth=0.8, label=label)

    ax.set_xlabel("time /s")
    ax.set_yscale("log")
    ax.legend(fontsize=13)

    ax.set_ylabel("droplet radius /\u03BCm")

    fig.tight_layout()

    if savename != "":
        fig.savefig(savename, dpi=400, bbox_inches="tight", facecolor="w", format="png")
        print("Figure .png saved as: " + str(savename))

    return fig, ax


def plot_randomsample_superdrops(time, superdrops, nsample, savename=""):
    """plot timeseries of the attributes of a
    random sample of superdroplets"""

    fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(10, 6))
    fig.suptitle("Time Evolution of a Random Sample of Superdroplets")

    superdrops.attach_time(time.mins, "min", do_reshape=True, var4reshape="sdId")

    sample_population = list(np.unique(ak.flatten(superdrops.sdId())))
    ids2plot = random.sample(sample_population, nsample)
    attrs = ["time", "radius", "xi", "msol", "coord3", "coord1", "coord2"]
    sample = superdrops.sample("sdId", sample_values=ids2plot, variables2sample=attrs)
    sample_time = sample.time()

    for a, attr in enumerate(["radius", "xi", "msol"]):
        data = sample[attr]
        if len(ak.flatten(data)) == 0:
            warnings.warn(f"no data for {attr} found. Not plotting {attr}")
            continue

        try:
            t = np.array(
                sample_time
            )  # will fail if time is different for different superdroplets
            d = np.array(data)
            axs[0, a].plot(t, d, linewidth=0.8)
        except ValueError:  # data not converible to numpy array (diff lengths of time for each superdroplet)?
            for i in range(len(data)):
                t = sample_time[i]
                d = data[i]
                axs[0, a].plot(t, d, linewidth=0.8)  # plot each SD seperately

    mks = MarkerStyle("o", fillstyle="full")
    for a, coord in enumerate(["coord3", "coord1", "coord2"]):
        data = sample[coord]
        if len(ak.flatten(data)) == 0:
            warnings.warn(f"no data for {coord} found. Not plotting {coord}")
            continue

        try:
            t = np.array(sample_time)
            d = np.array(data) / 1000  # [km]
            axs[1, a].plot(t, d, linestyle="", marker=mks, markersize=0.2)
        except (
            ValueError
        ):  # data not converible to numpy array (diff lengths of time for SDs)?
            for i in range(len(data)):  # plot each SD seperately
                t = sample_time[i]
                d = data[i] / 1000  # [km]
                axs[1, a].plot(t, d, linestyle="", marker=mks, markersize=0.2)

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

    return fig, axs


def plot_randomsample_superdrops_2dmotion(
    superdrops,
    nsample,
    savename="",
    arrows=False,
    ids2plot=None,
    israndom=True,
    cmap_var=[None, None],
    fig=None,
    ax=None,
):
    """plot timeseries of the attributes of a
    random sample of superdroplets"""
    if fig is None:
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6, 6))

    if ids2plot is None:
        if israndom:
            sample_population = list(np.unique(ak.flatten(superdrops.sdId())))
            ids2plot = random.sample(sample_population, nsample)
        else:
            minid = ak.min(superdrops.sdId())
            maxid = ak.max(superdrops.sdId())
            ids2plot = list(np.linspace(minid, maxid, nsample, dtype=int))

    mks = MarkerStyle("o", fillstyle="full")
    sample = superdrops.sample(
        "sdId", sample_values=ids2plot, variables2sample=["coord3", "coord1", "sdId"]
    )
    coordz = ak.flatten(sample.coord3()) / 1000  # [km]
    coordx = ak.flatten(sample.coord1()) / 1000  # [km]

    if all([c is None for c in cmap_var]):
        cmap = "gist_earth"
        c = ak.flatten(sample.sdId())
    else:
        cmap = cmap_var[0]
        c = superdrops.sample(
            "sdId", sample_values=ids2plot, variables2sample=cmap_var[1]
        )
        c = ak.flatten(c[cmap_var[1]])
    ax.scatter(coordx, coordz, marker=mks, s=0.4, cmap=cmap, c=c)

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

    return fig, ax
