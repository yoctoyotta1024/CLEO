"""
Copyright (c) 2025 MPI-M, Clara Bayley


----- CLEO -----
File: icon_bubble.py
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


# ICON
#
# ---------------------------------------------------------------
# Copyright (C) 2004-2025, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# Contact information: icon-model.org
#
# See AUTHORS.TXT for a list of authors
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ---------------------------------------------------------------

# %%
from functools import wraps
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.cm import ScalarMappable
import numpy as np
import pandas as pd
import seaborn as sns
import xarray as xr
import json
import datetime

MASS_CONSERVATION_THRESHOLD = 1e-10
DOMAIN_RESOLUTION_X = 5  # [km]
DOMAIN_SIZE = 100.0  # [km] (5km resolution => 40 triangles per latitude band)
DOMAIN_SIZE_Y = (DOMAIN_RESOLUTION_X * np.sqrt(3) / 2) * 4  # [km] (4 latitude bands)

from moist_thermodynamics.constants_icon import (
    Tmelt,
    ci,
    clw,
    cpd,
    cpv,
    g,
    ls,
    lv,
    rd,
    rv,
)

cvd = cpd - rd
cvv = cpv - rv
lf = ls - lv
lvc = lv - (cpv - clw) * Tmelt
lsc = ls - (cpv - ci) * Tmelt

plt.close("all")


# %%
def add_half_levels(ds):
    nlev = ds.sizes["height"]
    zh = xr.DataArray(
        np.zeros(nlev + 1),
        dims="height_2",
        attrs={"long_name": "altitude", "units": "km"},
    )
    for k in np.arange(nlev, 0, -1) - 1:
        zh[k] = 2 * ds.zg[k, 0] / 1000.0 - zh[k + 1]
    ds["zh"] = zh
    return ds.swap_dims({"height_2": "zh"})


def replace_icon2datetime(ds):
    def icon2datetime(icon_dates):
        """
        Note:
            - The function supports conversion of ICON format datetime values to Python datetime objects.
            - ICON format represents dates as numeric values (e.g., 20011201.5 for December 1, 2001, 12:00:00).
            - This function has been taken from the esm_analysis library: https://github.com/pgierz/esm_analysis
        """
        try:
            icon_dates = icon_dates.values
        except AttributeError:
            pass

        try:
            icon_dates = icon_dates[:]
        except TypeError:
            icon_dates = np.array([icon_dates])

        def _convert(icon_date):
            frac_day, date = np.modf(icon_date)
            frac_day *= 60**2 * 24
            date = str(int(date))
            date_str = datetime.datetime.strptime(date, "%Y%m%d")
            td = datetime.timedelta(seconds=int(frac_day.round(0)))
            return date_str + td

        conv = np.vectorize(_convert)
        try:
            out = conv(icon_dates)
        except TypeError:
            out = icon_dates
        if len(out) == 1:
            return pd.DatetimeIndex(out)[0]
        return pd.DatetimeIndex(out)

    return ds.assign_coords({"time": icon2datetime(ds.time)})


def read_bubble(expdir: Path, expname: str):
    ds2d = xr.open_dataset(
        expdir / f"{expname}_atm_2d_ml_20080801T000000Z.nc",
        engine="netcdf4",
    )
    ds3d = xr.open_dataset(
        expdir / f"{expname}_atm_3d_ml_20080801T000000Z.nc",
        engine="netcdf4",
    )

    return (
        ds2d.merge(ds3d)
        .rename_dims({"height_3": "zh"})
        .pipe(add_half_levels)
        .assign_coords(zg=lambda ds: ds.zg)
        .pipe(replace_icon2datetime)
    )


def get_crossection(ds):
    zg = (ds.zg[:, 0] / 1000).assign_attrs({"long_name": "altitude", "units": "km"})
    # conversion factor is from latitude in radians to km,
    # because 2*np.pi [radians] == (ds.clon_bnds.max().values - ds.clon_bnds.min().values) [radians] == DOMAIN_SIZE [km]
    conversion = DOMAIN_SIZE / (2 * np.pi)
    xg = xr.DataArray(
        ds.clon_bnds[:, 2] * conversion,
        dims="ncells",
        attrs={"long_name": "x", "units": "km"},
    )
    clat = xr.DataArray((ds.clat_bnds.mean(dim="vertices").values), dims="ncells")
    cond1 = clat == clat[0].values

    return (
        ds.assign_coords({"z": zg})
        .assign_coords({"x": xg})
        .where(cond1, drop=True)
        .swap_dims({"height": "z"})
        .swap_dims({"ncells": "x"})
        .sel(z=slice(6000.0, 0))
    )


def get_zy_crossection(ds):
    zg = (ds.zg[:, 0] / 1000).assign_attrs({"long_name": "altitude", "units": "km"})
    # conversion factor is from longitude in radians to km,
    # because 2*np.pi*5/160 [radians] == (ds.clat_bnds.max().values - ds.clat_bnds.min().values) [radians] == DOMAIN_SIZE_Y [km]
    conversion = DOMAIN_SIZE_Y / (2 * np.pi * 5 / 160)
    yg = xr.DataArray(
        ds.clat_bnds[:, 2] * conversion,
        dims="ncells",
        attrs={"long_name": "y", "units": "km"},
    )
    cell_idx = (
        ds.ncells.size // 2
    )  # x-section through approx. middle longitude of domain
    clon = xr.DataArray((ds.clon_bnds.mean(dim="vertices").values), dims="ncells")
    cond1 = clon == clon[cell_idx].values

    return (
        ds.assign_coords({"z": zg})
        .assign_coords({"y": yg})
        .where(cond1, drop=True)
        .swap_dims({"height": "z"})
        .swap_dims({"ncells": "y"})
        .sel(z=slice(6000.0, 0))
    )


def timestepped_plot(plot_function):
    @wraps(plot_function)
    def f(ds, t0, fig_axs=[], zyplane=False):
        if fig_axs == []:
            fig, axs = plt.subplots(
                1,
                len(ds.time),
                sharey=True,
                figsize=(12, 4),
                constrained_layout=True,
            )
        else:
            fig, axs = fig_axs
        for i, ax in enumerate(axs):
            sub_ds = ds.isel(time=i)

            plot_function(sub_ds, ax, zyplane)

            sns.despine(offset=0, left=True)
            ax.label_outer()
            ax.set_ylim(0, 6.1)
            if zyplane:
                ax.set_xlabel(f"{ds.y.long_name} / {ds.y.units}")
            else:
                ax.set_xlabel(f"{ds.x.long_name} / {ds.x.units}")
            time_minutes = (sub_ds.time - t0) / np.timedelta64(60, "s")
            ax.set_title(f"{time_minutes:.0f} min", fontsize=8)

        axs[0].set_ylabel(f"{ds.z.long_name} / {ds.z.units}")

        return fig

    return f


@timestepped_plot
def plot_condensate(ds, ax, zyplane):
    (ds.clw * 1000).plot.contourf(
        ax=ax,
        levels=[0.001, 0.4, 2.0, 4.0, 20.0],
        colors=["white", "lightgrey", "darkgrey", "grey", "black", "purple"],
        add_colorbar=False,
    )
    (ds.qr * 1000).plot.contour(
        ax=ax,
        levels=[0.001, 0.2, 1.5],
        colors=["darkblue", "darkblue", "darkblue"],
        linestyles=["dotted", "dashed", "solid"],
        linewidths=1.5,
    )
    (ds.qs * 10000).plot.contour(
        ax=ax,
        levels=[0.01, 0.1, 0.5],
        colors=["dodgerblue", "dodgerblue", "dodgerblue"],
        linestyles=["dotted", "dashed", "solid"],
        linewidths=1.5,
    )
    (ds.cli * 1000).plot.contourf(
        ax=ax,
        levels=[0.0, 0.01, 0.2, 1.5],
        colors=["white", "pink", "pink", "pink"],
        alpha=0.2,
        add_colorbar=False,
    )
    (ds.qg * 1000).plot.contour(
        ax=ax,
        levels=[0.01, 0.2, 1.5],
        colors=["fuchsia", "fuchsia", "fuchsia"],
        linestyles=["dotted", "dashed", "solid"],
        linewidths=1.5,
    )
    ax.set_xlim(-20.0, 20.0)


@timestepped_plot
def plot_temperature(ds, ax, zyplane):
    w = ds.wa
    qc = ds.clw + ds.cli + ds.qr + ds.qs + ds.qg
    Tv = ds.ta * (1 + 0.608 * ds.hus - qc)
    Tp = Tv - Tv.mean(dim="x")
    w.plot.contourf(ax=ax, cmap="bwr", levels=np.arange(-3, 3, 0.5), add_colorbar=False)
    Tp.plot.contour(
        ax=ax,
        levels=[-1.0, -0.5, 0.5, 1.0],
        colors=["dimgrey", "k"],
        linewidths=[1.0, 0.5, 0.5, 1.0],
        linestyles=["dashed", "dashed", "solid", "solid"],
    )
    ax.set_xlim(-50.0, 50.0)


def plot_winds(ds, t0, zyplane=False):
    @timestepped_plot
    def plot_wa(ds, ax, zyplane):
        if zyplane:
            w = ds.wa.sortby(ds.wa.y)
        else:
            w = ds.wa.sortby(ds.wa.x)
        w.plot.pcolormesh(ax=ax, cmap="bwr", vmin=-3, vmax=3, add_colorbar=False)

    @timestepped_plot
    def plot_ua(ds, ax, zyplane):
        if zyplane:
            u = ds.ua.sortby(ds.ua.y)
        else:
            u = ds.ua.sortby(ds.ua.x)
        u.plot.pcolormesh(ax=ax, cmap="bwr", vmin=-3, vmax=3, add_colorbar=False)

    @timestepped_plot
    def plot_va(ds, ax, zyplane):
        if zyplane:
            v = ds.va.sortby(ds.va.y)
        else:
            v = ds.va.sortby(ds.va.x)
        v.plot.pcolormesh(ax=ax, cmap="bwr", vmin=-3, vmax=3, add_colorbar=False)

    fig, axes = plt.subplots(
        3,
        len(ds.time) + 1,
        figsize=(16, 12),
        constrained_layout=True,
        width_ratios=[27] * len(ds.time) + [1],
    )
    axs = axes[:, :-1]
    caxs = axes[:, -1]

    plot_wa(ds, t0, fig_axs=[fig, axs[0, :]], zyplane=zyplane)
    plot_ua(ds, t0, fig_axs=[fig, axs[1, :]], zyplane=zyplane)
    plot_va(ds, t0, fig_axs=[fig, axs[2, :]], zyplane=zyplane)

    cmap = plt.get_cmap("bwr")
    norm = mcolors.Normalize(vmin=-3, vmax=3)
    for i, (cax, var) in enumerate(zip(caxs, ["wa", "ua", "va"])):
        axes[i, 0].set_title(f"{var}\n{axes[i,0].get_title()}")
        fig.colorbar(
            ScalarMappable(norm=norm, cmap=cmap),
            cax=cax,
            orientation="vertical",
            label=ds[var].attrs["long_name"],
        )

    for ax in axs.flatten():
        ax.set_xlim(-50.0, 50.0)
        ax.set_ylim(0, 3.0)

    fig.suptitle("ICON output of bubble winds")

    return fig


def plot_thermo(ds, t0, zyplane=False):
    @timestepped_plot
    def plot_press(ds, ax, zyplane):
        if zyplane:
            dp = (ds.pfull - ds.pfull.mean(dim="y")).sortby(ds.pfull.y)
        else:
            dp = (ds.pfull - ds.pfull.mean(dim="x")).sortby(ds.pfull.x)
        dp.plot.pcolormesh(ax=ax, cmap="plasma", vmin=-50, vmax=50, add_colorbar=False)

    @timestepped_plot
    def plot_temp(ds, ax, zyplane):
        if zyplane:
            ta = ds.ta.sortby(ds.ta.y)
        else:
            ta = ds.ta.sortby(ds.ta.x)
        ta.plot.pcolormesh(ax=ax, cmap="PuRd", vmin=270, vmax=290, add_colorbar=False)

    @timestepped_plot
    def plot_hus(ds, ax, zyplane):
        if zyplane:
            hus = ds.hus.sortby(ds.hus.y)
        else:
            hus = ds.hus.sortby(ds.hus.x)
        hus.plot.pcolormesh(
            ax=ax, cmap="PuBuGn", vmin=0.002, vmax=0.012, add_colorbar=False
        )

    fig, axes = plt.subplots(
        3,
        len(ds.time) + 1,
        figsize=(16, 12),
        constrained_layout=True,
        width_ratios=[27] * len(ds.time) + [1],
    )
    axs = axes[:, :-1]
    caxs = axes[:, -1]

    plot_press(ds, t0, fig_axs=[fig, axs[0, :]], zyplane=zyplane)
    plot_temp(ds, t0, fig_axs=[fig, axs[1, :]], zyplane=zyplane)
    plot_hus(ds, t0, fig_axs=[fig, axs[2, :]], zyplane=zyplane)

    # (!) must match plot_press, plot_temp and plot_hus vars, colormaps and vmin/vmax values with the colorbar
    vars = ["delta_press", "ta", "hus"]
    cmaps = ["plasma", "PuRd", "PuBuGn"]
    vmins = [-50, 270, 0.002]
    vmaxs = [50, 290, 0.012]

    for i, (cax, var, cmap, vmin, vmax) in enumerate(
        zip(caxs, vars, cmaps, vmins, vmaxs)
    ):
        axes[i, 0].set_title(f"{var}\n{axes[i,0].get_title()}")
        if var == "delta_press":
            label = f"$\u0394$ pressure / {ds['pfull'].attrs['units']}"
        else:
            label = f"{ds[var].attrs['long_name']} / {ds[var].attrs['units']}"
        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
        fig.colorbar(
            ScalarMappable(norm=norm, cmap=plt.get_cmap(cmap)),
            cax=cax,
            orientation="vertical",
            label=label,
        )

    for ax in axs.flatten():
        ax.set_xlim(-50.0, 50.0)
        ax.set_ylim(0, 3.0)

    return fig


def plot_liquid_water(ds, t0, zyplane=False):
    @timestepped_plot
    def plot_water(ds, ax, zyplane):
        norm = mcolors.LogNorm(vmin=0.001, vmax=10)
        if zyplane:
            qx = (ds.qr + ds.clw).sortby(ds.qr.y) * 1000
        else:
            qx = (ds.qr + ds.clw).sortby(ds.qr.x) * 1000
        qx.plot.pcolormesh(ax=ax, cmap="GnBu", norm=norm, add_colorbar=False)

    @timestepped_plot
    def plot_cloud(ds, ax, zyplane):
        norm = mcolors.LogNorm(vmin=0.001, vmax=10)
        if zyplane:
            qc = ds.clw.sortby(ds.clw.y) * 1000
        else:
            qc = ds.clw.sortby(ds.clw.x) * 1000
        qc.plot.pcolormesh(ax=ax, cmap="GnBu", norm=norm, add_colorbar=False)

    @timestepped_plot
    def plot_rain(ds, ax, zyplane):
        norm = mcolors.LogNorm(vmin=0.001, vmax=10)
        if zyplane:
            qr = ds.qr.sortby(ds.qr.y) * 1000
        else:
            qr = ds.qr.sortby(ds.qr.x) * 1000
        qr.plot.pcolormesh(ax=ax, cmap="GnBu", norm=norm, add_colorbar=False)

    fig, axs = plt.subplots(
        3,
        len(ds.time),
        sharey=True,
        figsize=(16, 12),
        constrained_layout=True,
    )
    plot_water(ds, t0, fig_axs=[fig, axs[0, :]], zyplane=zyplane)
    plot_cloud(ds, t0, fig_axs=[fig, axs[1, :]], zyplane=zyplane)
    plot_rain(ds, t0, fig_axs=[fig, axs[2, :]], zyplane=zyplane)
    for ax in axs.flatten():
        ax.set_xlim(-50.0, 50.0)
        ax.set_ylim(0, 3.0)

    return fig


def align_time(time, data):
    return xr.align(time, data, join="left")[1]


def attach_mass_summaries(ds):
    dz = xr.DataArray(-1000 * ds.zh.diff(dim="zh").values, dims="height")
    M = ds.rho * dz
    Mv = ds.hus * M
    Mc = ds.clw * M
    Mx = (ds.cli + ds.qr + ds.qg + ds.qs) * M
    Mt = Mv + Mc + Mx
    Md = M - Mt

    dM = align_time(ds.time, M.diff(dim="time").sum(axis=(1, 2))).fillna(0)
    dMt = align_time(ds.time, Mt.diff(dim="time").sum(axis=(1, 2))).fillna(0)
    Msfc = (ds.evspsbl + ds.prlr + ds.prls).sum(dim="ncells") * 30.0

    return ds.assign(
        dz=dz,
        M=M,
        Mv=Mv,
        Mc=Mc,
        Mx=Mx,
        Mt=Mt,
        Md=Md,
        dM=dM,
        dMt=dMt,
        Msfc=Msfc,
    )


def plot_mass(ds):
    ds = ds.pipe(attach_mass_summaries)

    sns.set_context("paper")
    fig, ax = plt.subplots(1, 1, sharey=True, figsize=(5, 3), constrained_layout=True)

    ds.dM.plot(label="$\\int \\rho \\, \\mathrm{d}z$", c="k")
    # ds.dMt.plot(label='$\\int q_\\mathrm{t} \\rho \\, \\mathrm{d}z$',c='k',ls='dashed')
    (ds.dMt + ds.Msfc).plot(
        label="$\\int q_\\mathrm{t} \\rho \\, \\mathrm{d}z + Q_\\mathrm{t} \\mathrm{d}t$",
        c="dodgerblue",
    )

    ax.set_ylabel("d$M$ / kgm$^{-2}$")
    ax.set_ylim(-MASS_CONSERVATION_THRESHOLD, MASS_CONSERVATION_THRESHOLD)
    plt.legend(ncol=1)
    sns.despine(offset=10)

    return fig


def plot_energy(ds):
    dt = ds.time.diff(dim="time") / np.timedelta64(1, "s")

    M = ds.rho * xr.DataArray(
        -ds.zh.diff(dim="zh").values, dims="height"
    )  # *1000.  # recall zh in km, indexed from top
    PE = (ds.zg[:, 0] * g * M).sum(dim="height").mean(dim="ncells")
    KE = (
        (0.5 * (ds.ua**2 + ds.va**2 + ds.wa.interp(zh=ds.zg.sel(ncells=0)) ** 2) * M)
        .sum(dim="height")
        .mean(dim="ncells")
    )

    Tk = ds.ta
    qv = ds.hus
    qliq = ds.clw + ds.qr
    qice = ds.cli + ds.qs + ds.qg
    qtot = qliq + qice + qv
    cx = cvd * (1 - qtot) + cvv * qv + clw * qliq + ci * qice
    IE = (
        (
            (
                cx * Tk
                - qliq * (lv - (cpv - clw) * Tmelt)
                - qice * (ls - (cpv - ci) * Tmelt)
            )
            * M
        )
        .mean(dim="ncells")
        .sum(dim="height")
    )

    dPEdt = PE.diff(dim="time") / dt
    dKEdt = KE.diff(dim="time") / dt
    dIEdt = IE.diff(dim="time") / dt
    dTEdt = dPEdt + dKEdt + dIEdt

    sns.set_context("paper")
    fig, ax = plt.subplots(1, 1, sharey=True, figsize=(5, 3), constrained_layout=True)

    dIEdt.plot(
        label="$\\mathrm{d}E_\\mathrm{I}/\\mathrm{d}t$",
        c="crimson",
        ls="dotted",
    )
    dPEdt.plot(
        label="$\\mathrm{d}E_\\mathrm{P}/\\mathrm{d}t$",
        c="dodgerblue",
        ls="dotted",
    )
    dKEdt.plot(label="$\\mathrm{d}E_\\mathrm{K}/\\mathrm{d}t$", c="green", ls="dotted")
    dTEdt.plot(label="$\\mathrm{d}E/\\mathrm{d}t$", c="k")

    ax.set_ylabel("Power density / Wm$^{-2}$")
    plt.legend(ncol=2)
    sns.despine(offset=10)

    return fig


# %%
def make_zx_plots(
    ds: xr.Dataset,
    ds_zx: xr.Dataset,
    outdir: Path,
    label: str = "",
    time_values: pd.DatetimeIndex = None,
):
    sub_ds = ds_zx.sel(time=time_values)

    plot_winds(sub_ds, ds.time[0]).savefig(outdir / f"winds{label}.png")
    plot_thermo(sub_ds, ds.time[0]).savefig(outdir / f"thermo{label}.png")
    plot_liquid_water(sub_ds, ds.time[0]).savefig(outdir / f"liquid_water{label}.png")
    plot_condensate(sub_ds, ds.time[0]).savefig(outdir / f"zzz_condensate{label}.png")
    plot_temperature(sub_ds, ds.time[0]).savefig(outdir / f"zzz_temperature{label}.png")
    plot_mass(ds).savefig(outdir / f"zzz_mass{label}.png")
    plot_energy(ds).savefig(outdir / f"zzz_energy{label}.png")


def plot_zy_winds(
    ds: xr.Dataset,
    ds_zy: xr.Dataset,
    outdir: Path,
    label: str = "",
    time_values: pd.DatetimeIndex = None,
):
    sub_ds = ds_zy.sel(time=time_values)

    plot_winds(sub_ds, ds.time[0], zyplane=True).savefig(
        outdir / f"winds_zyplane{label}.png"
    )
    plot_thermo(sub_ds, ds.time[0], zyplane=True).savefig(
        outdir / f"thermo_zyplane{label}.png"
    )
    plot_liquid_water(sub_ds, ds.time[0], zyplane=True).savefig(
        outdir / f"liquid_water_zyplane{label}.png"
    )

    return sub_ds


def get_metrics(ds):
    metrics = {}
    for varname in ["wa", "qs", "qr", "qg", "clw", "cli"]:
        metrics[varname] = {
            "min": float(ds[varname].min().values),
            "max": float(ds[varname].max().values),
        }
    return metrics


def write_metrics(ds: xr.Dataset, outdir: Path, label: str = ""):
    metrics = get_metrics(ds)
    with open(outdir / f"metrics{label}.json", "w") as outfile:
        json.dump(metrics, outfile, indent=4)


# %%
expname = "bubble_1mom"  # probably "bubble_1mom" or "bubble_cleo"
expdir = (
    Path("/work")
    / "mh0731"
    / "m300950"
    / "icon-mpim"
    / "experiments"
    / expname
    / "outdata"
)

outdir = expdir / "plots"
assert expdir.is_dir(), "expdir must be existing directory"
outdir.mkdir(exist_ok=True)

ds = read_bubble(expdir, expname)
print(f"EXPNAME: {expname}")
ds
# %%
write_metrics(ds, outdir)

# %% takes circa. 2mins
ds_zx = ds.pipe(get_crossection)
ds_zyplane = ds.pipe(get_zy_crossection)
# %%
time_values = pd.date_range(
    start="2008-08-01 00:40:00", end="2008-08-01 02:00:00", freq="20min"
)
make_zx_plots(ds, ds_zx, outdir, time_values=time_values)
plot_zy_winds(ds, ds_zyplane, outdir, time_values=time_values)
plt.close("all")
# %%
time_values = pd.date_range(
    start="2008-08-01 00:00:00", end="2008-08-01 00:40:00", freq="10min"
)
make_zx_plots(ds, ds_zx, outdir, label="_1", time_values=time_values)
plot_zy_winds(ds, ds_zyplane, outdir, label="_1", time_values=time_values)
plt.close("all")
# %%
print("longitudes: max, min, and: delta / radians = DOMAIN_SIZE / km")
print(
    f"max: {ds.clon_bnds.max().values}\nmin: {ds.clon_bnds.min().values}\n",
    f"2*pi*{(ds.clon_bnds.max().values - ds.clon_bnds.min().values) / (2 * np.pi)} radians = {DOMAIN_SIZE} km\n",
)
print("4 latitude bands")
print("latitudes: max, min, and: delta / radians = DOMAIN_SIZE_Y / km")
print(
    f"max: {ds.clat_bnds.max().values}\nmin: {ds.clat_bnds.min().values}\n",
    f"2*pi*{(ds.clat_bnds.max().values - ds.clat_bnds.min().values) / (2 * np.pi)} radians = {DOMAIN_SIZE_Y} km\n",
    f"(Hint: 5/160={5/160})",
)

# %%
time_values = pd.date_range(
    start="2008-08-01T01:58:00.000000000", end=ds.time.max().values, freq="30s"
)
make_zx_plots(ds, ds_zx, outdir, label="_2", time_values=time_values)
plot_zy_winds(ds, ds_zyplane, outdir, label="_2", time_values=time_values)
plt.close("all")
# %%
