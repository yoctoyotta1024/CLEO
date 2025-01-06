"""
Copyright (c) 2024 MPI-M, Clara Bayley


----- CLEO -----
File: bubble3d_plotting.py
Project: bubble3d
Created Date: Monday 6th January 2025
Author: Clara Bayley (CB)
Additional Contributors:
-----
Last Modified: Tuesday 7th January 2025
Modified By: CB
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
Script plots 3D example with time varying thermodynamics
for bubble test case output
"""

# %%
import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.cm import ScalarMappable
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument(
    "--path2CLEO",
    type=Path,
    help="Absolute path to CLEO",
    default="/home/m/m300950/CLEO",
)
parser.add_argument(
    "--path2build",
    type=Path,
    help="Absolute path to build",
    default="/home/m/m300950/CLEO/build_bubble3d",
)
parser.add_argument(
    "--config_filename",
    type=Path,
    help="Absolute path to config file",
    default="/home/m/m300950/CLEO/examples/bubble3d/src/config/bubble3d_config.yaml",
)
args, unknown = parser.parse_known_args()

sys.path.append(str(args.path2CLEO))  # imports from pySD
sys.path.append(
    str(args.path2CLEO / "examples" / "exampleplotting")
)  # imports from example plots package

from pySD.sdmout_src import pyzarr, pysetuptxt, pygbxsdat, sdtracing

# %%
path2build = args.path2build
config_filename = args.config_filename

dataset = path2build / "bin/bubble3d_sol.zarr"
setuptxt = path2build / "bin/bubble3d_setup.txt"
gridfile = path2build / "share/bubble3d_dimlessGBxboundaries.dat"
print(dataset)

ds = pyzarr.get_rawdataset(dataset)
config = pysetuptxt.get_config(setuptxt, nattrs=3, isprint=True)
consts = pysetuptxt.get_consts(setuptxt, isprint=True)
gbxs = pygbxsdat.get_gridboxes(gridfile, consts["COORD0"], isprint=True)

time = pyzarr.get_time(ds)
gbxindex = pyzarr.get_gbxindex(ds, gbxs["ndims"])
thermo, winds = pyzarr.get_thermodata(
    ds, config["ntime"], gbxs["ndims"], consts, getwinds=True
)
sddata = pyzarr.get_supers(ds, consts)
totnsupers = pyzarr.get_totnsupers(ds)

xfull_km = (gbxs["xfull"] - (gbxs["xfull"][-1] + gbxs["xfull"][0]) / 2) / 1000
zfull_km = gbxs["zfull"] / 1000
xxh_km, zzh_km = np.meshgrid(xfull_km, zfull_km, indexing="ij")


# %%
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

    savename = path2build / "bin" / f"bubble_{var}.png"
    fig.savefig(savename, dpi=400, bbox_inches="tight", facecolor="w", format="png")
    print("Figure .png saved as: " + str(savename))
# %%
maxid = np.amax(sddata.sdId[0])
nsample = 4799
sample_attrs = ["coord3", "coord1", "radius"]
sample = sdtracing.attrs_for_superdroplets_sample(
    sddata, sample_attrs, ndrops2sample=nsample, minid=0, maxid=maxid
)

# %%
wind_var = winds["wvel"]
label = "vertical velocity /m/s"
cmap = "PRGn"
vlims = [-5.0, 5.0]


# %%
def plot_2d_supers(xxh_km, zzh_km, wind_var, t2plts, sample, cmap, vlims):
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
    size = sample["radius"]
    coord3_km = sample["coord3"] / 1000
    coord1_km = sample["coord1"] / 1000 - 50

    for m in range(0, nplots):
        tidx = np.argmin(abs(time.mins - t2plts[m]))
        axs[m].contourf(xxh_km, zzh_km, wind_var[tidx, 0, :, :], cmap=cmap, norm=norm)
        axs[m].scatter(
            coord1_km[tidx, :], coord3_km[tidx, :], s=size[tidx, :], color="lightblue"
        )
        axs[m].set_title(
            "{:.0f} mins".format(time.mins[tidx]), fontsize=fontsize, y=1.025
        )
        axs[m].set_xlim([-50, 50])
        axs[m].set_ylim([0, 2.5])
        axs[m].spines["top"].set_visible(False)
        axs[m].spines["right"].set_visible(False)

    for ax in axs:
        ax.set_yticks([0, 1.25, 2.5])
        ax.set_yticklabels([0, "", 2.5], fontsize=fontsize)
        ax.set_xticks([-50, 0, 50])
        ax.set_xticklabels([-50, "", 50], fontsize=fontsize)
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


t2plts = [0, 30, 60, 90, 120]  # mins
fig, axs = plot_2d_supers(xxh_km, zzh_km, wind_var, t2plts, sample, cmap, vlims)
savename = path2build / "bin" / "bubble_motion.png"
fig.savefig(savename, dpi=400, bbox_inches="tight", facecolor="w", format="png")
print("Figure .png saved as: " + str(savename))
