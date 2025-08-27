"""
Copyright (c) 2024 MPI-M, Clara Bayley


----- CLEO -----
File: bubble3d_plotting.py
Project: bubble3d
Created Date: Monday 6th January 2025
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
Script plots 3D example with time varying thermodynamics
for bubble test case output
"""

# %%
### -------------------------------- IMPORTS ------------------------------- ###
import argparse
import awkward as ak
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.cm import ScalarMappable
from pathlib import Path


# %%
### ------------------------- FUNCTION DEFINITIONS ------------------------- ###
def parse_known_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--path2CLEO",
        type=Path,
        help="Absolute path to CLEO",
        default="/home/m/m300950/CLEO",
    )
    parser.add_argument(
        "--savefigpath",
        type=Path,
        help="Absolute path to build",
        default="/home/m/m300950/CLEO/build_bubble3d/bin",
    )
    parser.add_argument(
        "--grid_filename",
        type=Path,
        help="Absolute path to gridbox boundaries file",
        default="/home/m/m300950/CLEO/build_bubble3d/share/bubble3d_dimlessGBxboundaries.dat",
    )
    parser.add_argument(
        "--setupfile",
        type=Path,
        help="Absolute path to setup file",
        default="/home/m/m300950/CLEO/build_bubble3d/bin/bubble3d_setup.txt",
    )
    parser.add_argument(
        "--dataset",
        type=Path,
        help="Absolute path to dataset",
        default="/home/m/m300950/CLEO/build_bubble3d/bin/bubble3d_sol.zarr",
    )
    args, unknown = parser.parse_known_args()
    return args, unknown


def plot_2d_var(
    xxh_km, zzh_km, data, var, label, t2plts, cmap="plasma", vmin=None, vmax=None
):
    nplots = len(t2plts)
    fig, axes = plt.subplots(
        nrows=1,
        ncols=nplots + 1,
        figsize=(12, 4),
        constrained_layout=True,
        width_ratios=[27] * nplots + [1],
    )
    axs = axes[:-1]
    cax = axes[-1]

    cmap = plt.get_cmap(cmap)
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

    for m in range(0, nplots):
        tidx = np.argmin(abs(time.mins - t2plts[m]))
        axs[m].contourf(xxh_km, zzh_km, data[var][tidx, 0, :, :], cmap=cmap, norm=norm)
        axs[m].set_title("t={:.0f}mins".format(time.mins[tidx]), fontsize=10)
        axs[m].set_xlabel("x /km")
        axs[m].sharey(axs[0])
        axs[m].sharex(axs[0])
    axs[0].set_ylabel("z /km")

    fig.colorbar(
        ScalarMappable(norm=norm, cmap=cmap),
        cax=cax,
        orientation="vertical",
        label=label,
    )

    fig.suptitle("cleo_bubble: " + var)

    return fig, axes


def plot_2d_supers(xxh_km, zzh_km, wind_var, t2plts, sample, cmap, vlims, xlims):
    nplots = len(t2plts)
    fig, axes = plt.subplots(
        nrows=1,
        ncols=nplots + 1,
        figsize=(16, 5.5),
        constrained_layout=True,
        width_ratios=[24] * nplots + [1],
    )
    axs = axes[:-1]
    cax = axes[-1]
    fontsize = 18

    ## wind field contour plot
    cmap = plt.get_cmap(cmap)
    norm = mcolors.Normalize(vmin=vlims[0], vmax=vlims[1])

    ## superdroplets scatter plot
    sample_time = ak.Array(sample["time"])  # dims [superdrop, time(s)]
    radius = ak.Array(sample["radius"])
    coord3 = ak.Array(sample["coord3"])
    coord1 = ak.Array(sample["coord1"])
    shift = xlims[1]

    for m in range(0, nplots):
        t = t2plts[m]  # [mins]
        tidx = np.argmin(abs(time.mins - t))
        axs[m].contourf(xxh_km, zzh_km, wind_var[tidx, 0, :, :], cmap=cmap, norm=norm)

        delta = 1.0  # [mins]
        idx = ak.argmin(
            abs(t - sample_time), axis=1, keepdims=True
        )  # closest time(s) to t for each superdrop
        bad_time = (
            t - sample_time[idx] > delta
        )  # superdroplet times > delta away from time t
        size = ak.flatten(ak.where(bad_time, np.nan, radius[idx]))
        coord3_km = ak.flatten(ak.where(bad_time, np.nan, coord3[idx])) / 1000
        coord1_km = ak.flatten(ak.where(bad_time, np.nan, coord1[idx])) / 1000 - shift
        axs[m].scatter(coord1_km, coord3_km, s=size, color="lightblue")

        axs[m].set_title(
            "{:.0f} mins".format(time.mins[tidx]), fontsize=fontsize, y=1.025
        )
        axs[m].set_xlim(xlims)
        axs[m].set_ylim([0, 2.5])
        axs[m].spines["top"].set_visible(False)
        axs[m].spines["right"].set_visible(False)

    for ax in axs:
        ax.set_yticks([0, 1.25, 2.5])
        ax.set_yticklabels([0, "", 2.5], fontsize=fontsize)
        ax.set_xticks([xlims[0], 0, xlims[1]])
        ax.set_xticklabels([xlims[0], "", xlims[1]], fontsize=fontsize)
    for ax in axs[1:]:
        ax.set_xticklabels([])
        ax.set_yticklabels([])

    axs[0].set_xlabel("x /km", fontsize=fontsize)
    axs[0].set_ylabel("z /km", fontsize=fontsize)
    fig.suptitle(
        "Superdroplet motion in up/downdrafts modelled by CLEO-YAC-ICON",
        fontsize=fontsize,
        y=1.075,
    )
    cbar = fig.colorbar(
        ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="vertical"
    )
    cbar.set_label(label=label, fontsize=fontsize)
    ticks = [-5, 0, 5]
    cbar.set_ticks(ticks=ticks)
    cbar.set_ticklabels(ticklabels=ticks, fontsize=fontsize)

    return fig, axs


def plot_2d_supers_contours(
    xxh_km, zzh_km, wind_var, t2plts, sample, cmap, vlims, xlims
):
    nplots = len(t2plts)
    fig, axes = plt.subplots(
        nrows=1,
        ncols=nplots + 4,
        figsize=(18, 6),
        constrained_layout=True,
        width_ratios=[24] * nplots + [0.5, 1, 0.5, 1],
    )
    axs = axes[:-4]
    caxs = [axes[-3], axes[-1]]
    axes[-4].axis("off")
    axes[-2].axis("off")
    fontsize = 18

    ## wind field contour plot
    cmap = plt.get_cmap(cmap, 20)
    norm = mcolors.Normalize(vmin=vlims[0], vmax=vlims[1])

    ## superdroplets scatter plot
    sample_time = ak.Array(sample["time"])  # dims [superdrop, time(s)]
    radius = ak.Array(sample["radius"])
    coord3 = ak.Array(sample["coord3"])
    coord1 = ak.Array(sample["coord1"])
    shift = xlims[1]

    tidx = np.argmin(abs(time.mins - t2plts[-1]))
    sd_norm = mcolors.Normalize(vmin=time.mins[0], vmax=time.mins[tidx + 1])
    sd_cmap = "plasma"

    for m in range(0, nplots):
        t = t2plts[m]
        tidx = np.argmin(abs(time.mins - t))
        axs[m].contourf(xxh_km, zzh_km, wind_var[tidx, 0, :, :], cmap=cmap, norm=norm)

        idxs = sample_time <= t  # closest time(s) to t for each superdrop
        color = ak.flatten(ak.where(idxs, sample_time, np.nan))
        size = ak.flatten(ak.where(idxs, radius, np.nan)) * 5
        sd_y = ak.flatten(ak.where(idxs, coord3, np.nan)) / 1000
        sd_x = ak.flatten(ak.where(idxs, coord1, np.nan)) / 1000 - shift
        axs[m].scatter(
            sd_x,
            sd_y,
            s=size,
            c=color,
            cmap=sd_cmap,
            norm=sd_norm,
        )

        axs[m].set_title(
            "{:.0f} mins".format(time.mins[tidx]), fontsize=fontsize, y=1.025
        )
        axs[m].set_xlim(xlims)
        axs[m].set_ylim([0, 2.5])
        axs[m].spines["top"].set_visible(False)
        axs[m].spines["right"].set_visible(False)

    for ax in axs:
        ax.set_yticks([0, 1.25, 2.5])
        ax.set_yticklabels([0, "", 2.5], fontsize=fontsize)
        ax.set_xticks([xlims[0], 0, xlims[1]])
        ax.set_xticklabels([xlims[0], "", xlims[1]], fontsize=fontsize)
    for ax in axs[1:]:
        ax.set_xticklabels([])
        ax.set_yticklabels([])

    axs[0].set_xlabel("x /km", fontsize=fontsize, labelpad=-12)
    axs[0].set_ylabel("z /km", fontsize=fontsize, labelpad=-8)
    fig.suptitle(
        "Superdroplet motion in up/downdrafts modelled by CLEO-YAC-ICON",
        fontsize=fontsize + 2,
        y=1.075,
    )
    cbar = fig.colorbar(
        ScalarMappable(norm=norm, cmap=cmap), cax=caxs[0], orientation="vertical"
    )
    cbar.set_label(label=label, fontsize=fontsize, labelpad=-8)
    ticks = [-10, 0, 10]
    cbar.set_ticks(ticks=ticks)
    cbar.set_ticklabels(ticklabels=ticks, fontsize=fontsize)

    sd_cbar = fig.colorbar(
        ScalarMappable(norm=sd_norm, cmap=sd_cmap),
        cax=caxs[1],
        orientation="vertical",
    )
    sd_cbar.set_label(
        label="time of superdroplet location /mins", fontsize=fontsize, labelpad=-8
    )
    sd_ticks = [0, 20, 40, 60, 80]
    sd_ticklabels = [0, "", "", "", 80]
    sd_cbar.set_ticks(ticks=sd_ticks)
    sd_cbar.set_ticklabels(ticklabels=sd_ticklabels, fontsize=fontsize)

    return fig, axs


# %%
### --------------- IMPORT cleopy AND CLEO PLOTTING MODULES ------------------ ###
path2CLEO = parse_known_arguments()[0].path2CLEO
sys.path.append(
    str(path2CLEO / "examples" / "exampleplotting")
)  # imports from example plots package
from cleopy.sdmout_src import pyzarr, pysetuptxt, pygbxsdat

# %%
### -------------------------- INPUT PARAMETERS ---------------------------- ###
args = parse_known_arguments()[0]
savefigpath = args.savefigpath
dataset = args.dataset
setupfile = args.setupfile
grid_filename = args.grid_filename
print(dataset)

ds = pyzarr.get_rawdataset(dataset)
config = pysetuptxt.get_config(setupfile, nattrs=3, isprint=True)
consts = pysetuptxt.get_consts(setupfile, isprint=True)
gbxs = pygbxsdat.get_gridboxes(grid_filename, consts["COORD0"], isprint=True)

time = pyzarr.get_time(ds)
gbxindex = pyzarr.get_gbxindex(ds, gbxs["ndims"])
thermo, winds = pyzarr.get_thermodata(
    ds, config["ntime"], gbxs["ndims"], consts, getwinds=True
)
totnsupers = pyzarr.get_totnsupers(ds)
superdrops = pyzarr.get_supers(ds, consts)
superdrops.attach_time(time.mins, "min", do_reshape=True, var4reshape="sdId")

xfull_km = (gbxs["xfull"] - (gbxs["xfull"][-1] + gbxs["xfull"][0]) / 2) / 1000
zfull_km = gbxs["zfull"] / 1000
xxh_km, zzh_km = np.meshgrid(xfull_km, zfull_km, indexing="ij")

# %%
### -------------------------- CALL PLOT_2D_VAR ---------------------------- ###
t2plts = [0, 30, 60, 90, 120]  # mins
vars = ["wvel", "uvel", "vvel"]
labels = ["vertical velocity /m/s", "eastwards wind /m/s", "northwards wind /m/s"]
vlims = [[-3.0, 3.0], [-10.0, 10.0], [-3.0, 3.0]]
for v, var in enumerate(vars):
    vmin, vmax = vlims[v]
    label = labels[v]
    fig, axes = plot_2d_var(
        xxh_km, zzh_km, winds, var, label, t2plts, cmap="bwr", vmin=vmin, vmax=vmax
    )

    savename = savefigpath / f"bubble_{var}.png"
    fig.savefig(savename, dpi=400, bbox_inches="tight", facecolor="w", format="png")
    print("Figure .png saved as: " + str(savename))

# %%
### ------------ SAMPLE SDs AND SETTINGS FOR PLOT_2D_SUPERS ---------------- ###
nsample = 1000
sample_attrs = ["coord3", "coord1", "radius", "time"]
sample = superdrops.random_sample("sdId", nsample, variables2sample=sample_attrs)

wind_var = winds["wvel"]
label = "vertical velocity /m/s"
cmap = "PRGn"
vlims = [-5.0, 5.0]

# %%
### ------------------------ CALL PLOT_2D_SUPERS --------------------------- ###
t2plts = [0, 30, 60, 90, 120]  # mins
xl = (np.amax(gbxs["xhalf"]) - np.amin(gbxs["xhalf"])) / 2 / 1000
fig, axs = plot_2d_supers(
    xxh_km, zzh_km, wind_var, t2plts, sample, cmap, vlims, [-xl, xl]
)
savename = savefigpath / "bubble_motion.png"
fig.savefig(savename, dpi=400, bbox_inches="tight", facecolor="w", format="png")
print("Figure .png saved as: " + str(savename))

# %%
### -------- SAMPLE SDs AND SETTINGS FOR PLOT_2D_SUPERS_CONTOURS ----------- ###
nsample = 1000
sample_attrs = ["coord3", "coord1", "radius", "time"]
sample = superdrops.random_sample("sdId", nsample, variables2sample=sample_attrs)

wind_var = winds["wvel"]
label = "vertical velocity /m/s"
cmap = "PRGn"
vlims = [-10.0, 10.0]

# %%
### ------------------- CALL PLOT_2D_SUPERS_CONTOURS ----------------------- ###
t2plts = [20, 50, 80]  # mins
xl = (np.amax(gbxs["xhalf"]) - np.amin(gbxs["xhalf"])) / 2 / 1000
fig, axs = plot_2d_supers_contours(
    xxh_km, zzh_km, wind_var, t2plts, sample, cmap, vlims, [-xl, xl]
)
savename = savefigpath / "bubble_motion_v2.png"
fig.savefig(savename, dpi=400, bbox_inches="tight", facecolor="w", format="png")
print("Figure .png saved as: " + str(savename))
