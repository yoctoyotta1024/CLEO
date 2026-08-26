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
import matplotlib.pyplot as plt
import examples.bubble3d.bubble3d_plotting as bubbleplot
from importlib import reload
from cleopy.sdmout_src import pyzarr

reload(bubbleplot)

plt.close("all")
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

xhalf_km = (gbxs["xhalf"] - (gbxs["xhalf"][-1] + gbxs["xhalf"][0]) / 2) / 1000
zhalf_km = gbxs["zhalf"] / 1000
xxh_km, zzh_km = np.meshgrid(xhalf_km, zhalf_km, indexing="ij")

yhalf_km = (gbxs["yhalf"] - (gbxs["yhalf"][-1] + gbxs["yhalf"][0]) / 2) / 1000
yyh_km, zzyh_km = np.meshgrid(yhalf_km, zhalf_km, indexing="ij")

massmoms = pyzarr.get_massmoms(dataset, config["ntime"], gbxs["ndims"])
rainmassmoms = pyzarr.get_rainmassmoms(dataset, config["ntime"], gbxs["ndims"])

# %%
### -------------------------- CALL PLOT_2D_WINDS ---------------------------- ###
t2plts_0 = [40, 60, 80, 100, 120]  # mins
t2plts_1 = [0, 10, 20, 30, 40]  # mins
t2plts_2 = [118, 118.5, 119, 119.5, 120]  # mins

# %%
vars = ["wvel", "uvel", "vvel"]
labels = ["vertical velocity /m/s", "eastwards wind /m/s", "northwards wind /m/s"]
vlims = [[-3.0, 3.0]] * 3
for t2plts, figlabel in zip([t2plts_0, t2plts_1, t2plts_2], ["", "_1", "_2"]):
    nplots = len(t2plts)
    fig, axes = plt.subplots(
        nrows=3,
        ncols=nplots + 1,
        figsize=(16, 12),
        constrained_layout=True,
        width_ratios=[27] * nplots + [1],
    )
    for v, var in enumerate(vars):
        axs = axes[v, :-1]
        cax = axes[v, -1]
        vmin, vmax = vlims[v]
        label = labels[v]
        fig, axs = bubbleplot.plot_2d_winds(
            fig,
            axs,
            cax,
            xxh_km,
            zzh_km,
            winds,
            var,
            label,
            t2plts,
            cmap="bwr",
            vmin=vmin,
            vmax=vmax,
        )

    savename = savefigpath / f"bubble_winds{figlabel}.png"
    fig.savefig(savename, dpi=400, bbox_inches="tight", facecolor="w", format="png")
    print("Figure .png saved as: " + str(savename))
plt.close("all")
# %%
### ------------------- WIND PLOTS IN Z-Y (MID)PLANE ----------------------- ###
for t2plts, figlabel in zip([t2plts_0, t2plts_1], ["", "_1"]):
    nplots = len(t2plts)
    fig, axes = plt.subplots(
        nrows=3,
        ncols=nplots + 1,
        figsize=(16, 12),
        constrained_layout=True,
        width_ratios=[27] * nplots + [1],
    )

    for v, var in enumerate(vars):
        axs = axes[v, :-1]
        cax = axes[v, -1]
        vmin, vmax = vlims[v]
        label = labels[v]
        fig, axs = bubbleplot.plot_2d_var_zymidplane(
            fig,
            axs,
            cax,
            yyh_km,
            zzyh_km,
            winds,
            var,
            label,
            t2plts,
            cmap="bwr",
            vmin=vmin,
            vmax=vmax,
        )

        savename = savefigpath / f"bubble_zyplane_winds{figlabel}.png"
    fig.savefig(savename, dpi=400, bbox_inches="tight", facecolor="w", format="png")
    print("Figure .png saved as: " + str(savename))
plt.close("all")
# %%
### ----------------------------- PLOT THERMO ------------------------------ ###
dp = (
    thermo["press"][:, :, :, :] - np.mean(thermo["press"], axis=(2))[:, :, None, :]
) * 100
specific_humidity = thermo["qvap"] / 1000 / (1 + thermo["qvap"] / 1000)
data = {
    "delta_press": dp,
    "temp": thermo["temp"],
    "specific_humidity": specific_humidity,
}
labels = ["$\u0394$ pressure /Pa", "temperature /K", "specific humidity kg/kg"]
cmaps = ["plasma", "PuRd", "PuBuGn"]
vlims = [[-50, 50], [270, 290], [0.002, 0.012]]

for t2plts, figlabel in zip([t2plts_0, t2plts_1, t2plts_2], ["", "_1", "_2"]):
    nplots = len(t2plts)
    fig, axes = plt.subplots(
        nrows=3,
        ncols=nplots + 1,
        figsize=(16, 12),
        constrained_layout=True,
        width_ratios=[27] * nplots + [1],
    )
    for v, var in enumerate(list(data.keys())):
        axs = axes[v, :-1]
        cax = axes[v, -1]
        vmin, vmax = vlims[v]
        label = labels[v]
        fig, axs = bubbleplot.plot_2d_var(
            fig,
            axs,
            cax,
            xxh_km[:-1, :-1],
            zzh_km[:-1, :-1],
            data,
            var,
            label,
            t2plts,
            cmap=cmaps[v],
            vmin=vmin,
            vmax=vmax,
        )

    savename = savefigpath / f"bubble_thermo{figlabel}.png"
    fig.savefig(savename, dpi=400, bbox_inches="tight", facecolor="w", format="png")
    print("Figure .png saved as: " + str(savename))
plt.close("all")
# %%
### ------------------- WIND PLOTS IN Z-Y (MID)PLANE ----------------------- ###
for t2plts, figlabel in zip([t2plts_0, t2plts_1], ["", "_1"]):
    nplots = len(t2plts)
    fig, axes = plt.subplots(
        nrows=3,
        ncols=nplots + 1,
        figsize=(16, 12),
        constrained_layout=True,
        width_ratios=[27] * nplots + [1],
    )

    for v, var in enumerate(list(data.keys())):
        axs = axes[v, :-1]
        cax = axes[v, -1]
        vmin, vmax = vlims[v]
        label = labels[v]
        fig, axs = bubbleplot.plot_2d_var_zymidplane(
            fig,
            axs,
            cax,
            yyh_km,
            zzyh_km,
            data,
            var,
            label,
            t2plts,
            cmap=cmaps[v],
            vmin=vmin,
            vmax=vmax,
        )

        savename = savefigpath / f"bubble_zyplane_thermo{figlabel}.png"
    fig.savefig(savename, dpi=400, bbox_inches="tight", facecolor="w", format="png")
    print("Figure .png saved as: " + str(savename))
plt.close("all")
# %%
### ---------------------------- PLOT MASSMOMS ----------------------------- ###
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
for t2plts, figlabel in zip([t2plts_0, t2plts_1], ["", "_1"]):
    nplots = len(t2plts)
    fig, axes = plt.subplots(
        nrows=4,
        ncols=nplots + 1,
        figsize=(16, 16),
        constrained_layout=True,
        width_ratios=[27] * nplots + [1],
    )
    for v, var in enumerate(list(data.keys())):
        axs = axes[v, :-1]
        cax = axes[v, -1]
        vmin, vmax = vlims[v]
        label = labels[v]
        fig, axs = bubbleplot.plot_2d_var(
            fig,
            axs,
            cax,
            xxh_km[:-1, :-1],
            zzh_km[:-1, :-1],
            data,
            var,
            label,
            t2plts,
            cmap=cmaps[v],
            vmin=vmin,
            vmax=vmax,
        )

    savename = savefigpath / f"bubble_massmoms{figlabel}.png"
    fig.savefig(savename, dpi=400, bbox_inches="tight", facecolor="w", format="png")
    print("Figure .png saved as: " + str(savename))
plt.close("all")
# %%
### ------------ SAMPLE SDs AND SETTINGS FOR PLOT_2D_SUPERS ---------------- ###
nsample = 1440
sample_attrs = ["coord3", "coord1", "radius", "time"]
sample = superdrops.random_sample("sdId", nsample, variables2sample=sample_attrs)

wind_var = winds["wvel"]
label = "vertical velocity /m/s"
cmap = "PRGn"
vlims = [-3.0, 3.0]

# %%
### ------------------------ CALL PLOT_2D_SUPERS --------------------------- ###
t2plts = [30, 40, 50, 60, 70, 80, 90, 100, 110, 120]  # mins
xl = (np.amax(gbxs["xhalf"]) - np.amin(gbxs["xhalf"])) / 2 / 1000
fig, axs = bubbleplot.plot_2d_supers(
    xxh_km, zzh_km, wind_var, label, t2plts, sample, cmap, vlims, [-xl, xl]
)
savename = savefigpath / "bubble_motion.png"
fig.savefig(savename, dpi=400, bbox_inches="tight", facecolor="w", format="png")
print("Figure .png saved as: " + str(savename))
plt.close("all")
# %%
