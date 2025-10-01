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
