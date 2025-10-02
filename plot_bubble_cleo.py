"""
Copyright (c) 2025 MPI-M, Clara Bayley


----- CLEO -----
File: plot_bubble_cleo.py
Project: CLEO
Created Date: Thursday 1st January 1970
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
"""

# %%
### --------------------- IMPORTS AND DEFAULT PLOTS ------------------------ ###
import numpy as np
import examples.bubble3d.bubble3d_plotting as bubbleplot
from importlib import reload
from cleopy.sdmout_src import pyzarr

reload(bubbleplot)
# %%
print("zfull in CLEO /km:\n", bubbleplot.zfull_km)

# %%
### -------------------------- INPUT PARAMETERS ---------------------------- ###
savefigpath = bubbleplot.savefigpath
dataset = bubbleplot.dataset
setupfile = bubbleplot.setupfile
grid_filename = bubbleplot.grid_filename
print(dataset)

ds = bubbleplot.ds
config = bubbleplot.config
consts = bubbleplot.consts
gbxs = bubbleplot.gbxs

time = bubbleplot.time
gbxindex = bubbleplot.gbxindex
thermo, winds = bubbleplot.thermo, bubbleplot.winds
totnsupers = bubbleplot.totnsupers
superdrops = bubbleplot.superdrops
superdrops.attach_time(time.mins, "min", do_reshape=True, var4reshape="sdId")

xfull_km = bubbleplot.xfull_km
zfull_km = bubbleplot.zfull_km
xxh_km, zzh_km = bubbleplot.xxh_km, bubbleplot.zzh_km

yfull_km = (gbxs["yfull"] - (gbxs["yfull"][-1] + gbxs["yfull"][0]) / 2) / 1000
yyh_km, zzyh_km = np.meshgrid(yfull_km, zfull_km, indexing="ij")

massmoms = pyzarr.get_massmoms(dataset, config["ntime"], gbxs["ndims"])
rainmassmoms = pyzarr.get_rainmassmoms(dataset, config["ntime"], gbxs["ndims"])

# %%
### ----------------------------- PLOT THERMO ------------------------------ ###
t2plts = [30, 40, 50, 60, 70, 80, 90, 100, 110, 120]  # mins
dp = (
    thermo["press"][:, :, :, :] - np.mean(thermo["press"], axis=(2))[:, :, None, :]
) * 100
specific_humidity = thermo["qvap"] / 1000 / (1 + thermo["qvap"] / 1000)
data = {
    "press": dp,
    "temp": thermo["temp"],
    "specific_humidity": specific_humidity,
}
labels = ["pressure /Pa", "temperature /K", "specific humidity kg/kg"]
cmaps = ["plasma", "PuRd", "PuBuGn"]
vlims = [[-50, 50], [270, 290], [0.002, 0.012]]
for v, var in enumerate(list(data.keys())):
    vmin, vmax = vlims[v]
    label = labels[v]
    fig, axes = bubbleplot.plot_2d_var(
        xxh_km, zzh_km, data, var, label, t2plts, cmap=cmaps[v], vmin=vmin, vmax=vmax
    )

    savename = savefigpath / f"bubble_{var}.png"
    fig.savefig(savename, dpi=400, bbox_inches="tight", facecolor="w", format="png")
    print("Figure .png saved as: " + str(savename))

# %%
### ---------------------------- PLOT MASSMOMS ----------------------------- ###
t2plts = [30, 40, 50, 60, 70, 80, 90, 100, 110, 120]  # mins
data = {
    "massmom0": massmoms.mom0,
    "massmom1": massmoms.mom1,
    "rainmassmom0": rainmassmoms.mom0,
    "rainmassmom1": rainmassmoms.mom1,
}
labels = [
    "Number of Droplets",
    "Mass of Droplets /g",
    "Number of Raindrops",
    "Mass of Raindrops \g",
]
cmaps = ["GnBu"] * 4
vlims = [[0.0, 1e18], [0.0, 1e7]] * 2
for v, var in enumerate(list(data.keys())):
    vmin, vmax = vlims[v]
    label = labels[v]
    fig, axes = bubbleplot.plot_2d_var(
        xxh_km, zzh_km, data, var, label, t2plts, cmap=cmaps[v], vmin=vmin, vmax=vmax
    )

    savename = savefigpath / f"bubble_{var}.png"
    fig.savefig(savename, dpi=400, bbox_inches="tight", facecolor="w", format="png")
    print("Figure .png saved as: " + str(savename))

# %%
### ------------------ FUNCTION TO PLOT IN Z-Y (MID)PLANE ------------------ ###
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.cm import ScalarMappable


def plot_2d_var_zymidplane(
    yyh_km, zzh_km, data, var, label, t2plts, cmap="plasma", vmin=None, vmax=None
):
    nplots = len(t2plts)
    fig, axes = plt.subplots(
        nrows=1,
        ncols=nplots + 1,
        figsize=(16, 4),
        constrained_layout=True,
        width_ratios=[27] * nplots + [1],
    )
    axs = axes[:-1]
    cax = axes[-1]

    cmap = plt.get_cmap(cmap)
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

    xidx = np.shape(data[var])[2] // 2  # (middle coord or domain)
    if np.shape(data[var])[2] / 2 != xidx:
        print("WARNING: maybe not plotting middle of domain")

    for m in range(0, nplots):
        tidx = np.argmin(abs(time.mins - t2plts[m]))
        axs[m].pcolormesh(
            yyh_km, zzh_km, data[var][tidx, :, xidx, :], cmap=cmap, norm=norm
        )
        axs[m].set_title("t={:.0f}mins".format(time.mins[tidx]), fontsize=10)
        axs[m].set_xlabel("y /km")
        axs[m].sharey(axs[0])
        axs[m].sharex(axs[0])
        axs[m].spines["top"].set_visible(False)
        axs[m].spines["right"].set_visible(False)
        axs[m].spines["left"].set_visible(False)
    axs[0].set_ylabel("z /km")
    axs[0].set_ylim([0, 3.0])

    fig.colorbar(
        ScalarMappable(norm=norm, cmap=cmap),
        cax=cax,
        orientation="vertical",
        label=label,
    )

    fig.suptitle("cleo_bubble: " + var)

    return fig, axes


# %%
### ------------------- WIND PLOTS IN Z-Y (MID)PLANE ----------------------- ###
t2plts = [30, 40, 50, 60, 70, 80, 90, 100, 110, 120]  # mins
vars = ["wvel", "uvel", "vvel"]
labels = ["vertical velocity /m/s", "eastwards wind /m/s", "northwards wind /m/s"]
vlims = [[-3.0, 3.0]] * 3
for v, var in enumerate(vars):
    vmin, vmax = vlims[v]
    label = labels[v]
    fig, axes = plot_2d_var_zymidplane(
        yyh_km, zzyh_km, winds, var, label, t2plts, cmap="bwr", vmin=vmin, vmax=vmax
    )

    savename = savefigpath / f"bubble_zyplane_{var}.png"
    fig.savefig(savename, dpi=400, bbox_inches="tight", facecolor="w", format="png")
    print("Figure .png saved as: " + str(savename))

# %%
