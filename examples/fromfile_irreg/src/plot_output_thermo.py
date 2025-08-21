"""
Copyright (c) 2024 MPI-M, Clara Bayley


----- CLEO -----
File: plot_output_thermo.py
Project: src
Created Date: Wednesday 11th September 2024
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
Python functions used to make plots of CLEO's thermodynamics output for
fromfile_irreg example.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib import colors
from pathlib import Path


def plot_domain_thermodynamics_timeseries(time, gbxs, thermo, winds, savedir, do_show):
    """plot 2-D cross-sections of domain along y axis for thermodynamics and
    wind fields at a series of time slices."""

    labels = {
        "press": "Pressure /hPa",
        "temp": "Temperature /K",
        "qvap": "q$_{v}$ g/Kg",
        "qcond": "q$_{c}$ g/Kg",
        "wvel": "w m/s",
        "uvel": "u m/s",
        "vvel": "v m/s",
    }

    cmaps = {
        "press": ["PRGn", colors.CenteredNorm(vcenter=1015)],
        "temp": ["RdBu_r", colors.CenteredNorm(vcenter=300)],
        "qvap": ["BrBG", colors.CenteredNorm(vcenter=0.01)],
        "qcond": ["BrBG", colors.CenteredNorm(vcenter=0.001)],
        "wvel": ["coolwarm", colors.CenteredNorm(vcenter=0.0, halfrange=3.0)],
        "uvel": ["coolwarm", colors.CenteredNorm(vcenter=0.0, halfrange=3.0)],
        "vvel": ["coolwarm", colors.CenteredNorm(vcenter=0.0, halfrange=3.0)],
    }

    xxh, zzh = np.meshgrid(
        gbxs["xhalf"], gbxs["zhalf"], indexing="ij"
    )  # dims [xdims, zdims]
    t2plts = np.linspace(time.mins[0], time.mins[-1], 5)

    for key in labels.keys():
        try:
            data4d = thermo[key]
        except ValueError:
            data4d = winds[key]

        plot_timeseries_domain_slices(
            key,
            labels[key],
            cmaps[key][0],
            cmaps[key][1],
            time,
            t2plts,
            xxh,
            zzh,
            gbxs["yfull"],
            data4d,
            savedir,
        )
        if do_show:
            plt.show()


def plot_timeseries_domain_slices(
    key, label, cmap, norm, time, t2plts, xxh, zzh, yfull, data4d, savedir
):
    """plot 2-D cross-sections along y axis of 3-D data 'data3d' for several timeslices
    given times and meshgrid for centres in x and z, 'xxh' and 'zzh'."""

    fig = plt.figure(figsize=(16, 21))
    ncols = len(t2plts)
    nrows = data4d.shape[1]
    gs = GridSpec(
        nrows=nrows + 1, ncols=ncols, figure=fig, height_ratios=[1] + [8] * nrows
    )
    cax = fig.add_subplot(gs[0, :])

    cax.set_title(label)
    pcm = plot_domain_slices(
        fig, gs, nrows, ncols, time, t2plts, xxh, zzh, yfull, data4d, cmap, norm
    )
    plt.colorbar(pcm, cax=cax, orientation="horizontal")

    fig.tight_layout()

    if savedir != "":
        savename = savedir / Path(f"timeseries_{key}.png")
        fig.savefig(savename, dpi=400, bbox_inches="tight", facecolor="w", format="png")
        print("Figure .png saved as: " + str(savename))


def plot_domain_slices(
    fig, gs, nrows, ncols, time, t2plts, xxh, zzh, yfull, data4d, cmap, norm
):
    """plot 2-D cross-sections along y axis of 4-D data 'data4d' at a series of time slices
    given meshgrid for centres in x and z, 'xxh' and 'zzh'."""

    for m in range(ncols):
        t = np.argmin(abs(time.mins - t2plts[m]))
        data3d = data4d[t, :, :, :]
        for n in range(nrows):
            ax = fig.add_subplot(gs[n + 1, m])
            ax.set_title(
                "t={:.0f}min,\ny={:.2f}km".format(time.mins[t], yfull[n] / 1000)
            )
            pcm = ax_colormap(ax, xxh, zzh, data3d[n, :, :], cmap, norm)

    return pcm


def ax_colormap(ax, xxh, zzh, data2d, cmap, norm):
    """plot pcolormesh for 2-D data 'data2d' given meshgrid for
    centres in x and z, 'xxh' and 'zzh'."""

    ax.set_aspect("equal")
    pcm = ax.pcolormesh(
        xxh[:, :] / 1000, zzh[:, :] / 1000, data2d, cmap=cmap, norm=norm
    )
    ax.set_xlabel("x /km")
    ax.set_ylabel("z /km")

    return pcm
