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
import numpy as np
import pandas as pd
import seaborn as sns
import xarray as xr
import json
import datetime

MASS_CONSERVATION_THRESHOLD = 1e-10

# load thermodynamic variables

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
        .pipe(add_half_levels)
        .assign_coords(zg=lambda ds: ds.zg)
        .pipe(replace_icon2datetime)
    )


def get_crossection(ds):
    domain_size = 100.0
    zg = (ds.zg[:, 0] / 1000).assign_attrs({"long_name": "altitude", "units": "km"})
    xg = xr.DataArray(
        ds.clon_bnds[:, 2] / np.pi / 2.0 * domain_size,
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


def timestepped_plot(plot_function):
    @wraps(plot_function)
    def f(ds, t0, fig_axs=[]):
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

            plot_function(sub_ds, ax)

            sns.despine(offset=0, left=True)
            ax.label_outer()
            ax.set_ylim(0, 6.1)
            ax.set_xlabel(f"{ds.x.long_name} / {ds.x.units}")
            time_minutes = (sub_ds.time - t0) / np.timedelta64(60, "s")
            ax.set_title(f"{time_minutes:.0f} min", fontsize=8)

        axs[0].set_ylabel(f"{ds.z.long_name} / {ds.x.units}")

        return fig

    return f


@timestepped_plot
def plot_condensate(ds, ax):
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
def plot_temperature(ds, ax):
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


def plot_winds(ds, t0):
    @timestepped_plot
    def plot_wa(ds, ax):
        ds.wa.plot.contourf(
            ax=ax, cmap="bwr", levels=np.arange(-3, 3, 0.5), add_colorbar=False
        )

    @timestepped_plot
    def plot_ua(ds, ax):
        ds.ua.plot.contourf(
            ax=ax, cmap="bwr", levels=np.arange(-3, 3, 0.5), add_colorbar=False
        )

    @timestepped_plot
    def plot_va(ds, ax):
        ds.va.plot.contourf(
            ax=ax, cmap="bwr", levels=np.arange(-3, 3, 0.5), add_colorbar=False
        )

    fig, axs = plt.subplots(
        3,
        len(ds.time),
        sharey=True,
        figsize=(12, 8),
        constrained_layout=True,
    )
    plot_wa(ds, t0, fig_axs=[fig, axs[0, :]])
    plot_ua(ds, t0, fig_axs=[fig, axs[1, :]])
    plot_va(ds, t0, fig_axs=[fig, axs[2, :]])
    for ax in axs.flatten():
        ax.set_xlim(-50.0, 50.0)

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
def make_plots(ds: xr.Dataset, outdir: Path):
    time_values = pd.date_range(
        start="2008-08-01 00:30:00", end="2008-08-01 02:00:00", freq="10min"
    )
    sub_ds = ds.pipe(get_crossection).sel(time=time_values)

    plot_condensate(sub_ds, ds.time[0]).savefig(outdir / "condensate.png")
    plot_temperature(sub_ds, ds.time[0]).savefig(outdir / "temperature.png")
    plot_winds(sub_ds, ds.time[0]).savefig(outdir / "winds.png")
    plot_mass(ds).savefig(outdir / "mass.png")
    plot_energy(ds).savefig(outdir / "energy.png")


def get_metrics(ds):
    metrics = {}
    for varname in ["wa", "qs", "qr", "qg", "clw", "cli"]:
        metrics[varname] = {
            "min": float(ds[varname].min().values),
            "max": float(ds[varname].max().values),
        }
    return metrics


def write_metrics(ds, outdir):
    metrics = get_metrics(ds)
    with open(outdir / "metrics.json", "w") as outfile:
        json.dump(metrics, outfile, indent=4)


# %%
expname = "aes_bubble_cleo"
expdir = Path("/work/bm1183/m300950/icon-mpim/build/experiments/") / expname
outdir = expdir / "plots"
assert expdir.is_dir(), "expdir must be existing directory"
outdir.mkdir(exist_ok=True)

ds = read_bubble(expdir, expname)
# %%
make_plots(ds, outdir)
write_metrics(ds, outdir)
# %%
