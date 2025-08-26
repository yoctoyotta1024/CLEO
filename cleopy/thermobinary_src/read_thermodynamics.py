"""
Copyright (c) 2024 MPI-M, Clara Bayley


----- CLEO -----
File: read_thermodynamics.py
Project: thermobinary_src
Created Date: Friday 13th October 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from pathlib import Path

from .create_thermodynamics import thermoinputsdict, DimlessThermodynamics
from .thermogen import saturation_press
from ..readbinary import readbinary
from ..gbxboundariesbinary_src import read_gbxboundaries as rgrid


def relative_humidity(press, temp, qvap, Mr_ratio):
    pv = qvap * press / (Mr_ratio + qvap)  # vapour pressure
    psat = saturation_press(temp)
    relh = pv / psat

    qsat = Mr_ratio * psat / (press - pv)
    supersat = qvap / qsat - 1

    return relh, supersat


def potential_temperature(press, temp, press_ref, RGAS, CP):
    theta = temp * (press_ref / press) ** (RGAS / CP)

    return theta


class ThermoOnGrid:
    def __init__(self, thermodata, inputs, ndims):
        dth = DimlessThermodynamics(inputs=inputs)
        thermodata = dth.redimensionalise(thermodata)

        ntime = inputs["ntime"]
        ndims = [int(n) for n in ndims]
        cen = [inputs["ntime"]] + list(np.flip(ndims))
        zface = [ntime, ndims[2], ndims[1], ndims[0] + 1]
        xface = [ntime, ndims[2], ndims[1] + 1, ndims[0]]
        yface = [ntime, ndims[2] + 1, ndims[1], ndims[0]]

        self.vars = ["press", "temp", "qvap", "qcond"]
        self.press = np.reshape(thermodata["press"], cen)
        self.temp = np.reshape(thermodata["temp"], cen)
        self.qvap = np.reshape(thermodata["qvap"], cen)
        self.qcond = np.reshape(thermodata["qcond"], cen)

        if "wvel" in thermodata.keys():
            self.wvel = np.reshape(thermodata["wvel"], zface)
            self.wvel_cens = (self.wvel[:, :, :, 1:] + self.wvel[:, :, :, :-1]) / 2
            self.vars += ["wvel", "wvel_cens"]
            if "uvel" in thermodata.keys():
                self.uvel = np.reshape(thermodata["uvel"], xface)
                self.uvel_cens = (self.uvel[:, :, 1:, :] + self.uvel[:, :, :-1, :]) / 2
                self.vars += ["uvel", "uvel_cens"]
                if "vvel" in thermodata.keys():
                    self.vvel = np.reshape(thermodata["vvel"], yface)
                    self.vvel_cens = (
                        self.vvel[:, 1:, :, :] + self.vvel[:, :-1, :, :]
                    ) / 2
                    self.vars += ["vvel", "vvel_cens"]

    def __getitem__(self, key):
        if key not in self.vars:
            err = "no variable " + key + " in ThermoOnGrid"
            raise ValueError(err)
        elif key == "press":
            return self.press
        elif key == "temp":
            return self.temp
        elif key == "qvap":
            return self.qvap
        elif key == "qcond":
            return self.qcond

        elif key == "wvel":
            return self.wvel
        elif key == "uvel":
            return self.uvel
        elif key == "vvel":
            return self.vvel

        elif key == "wvel_cens":
            return self.wvel_cens
        elif key == "uvel_cens":
            return self.uvel_cens
        elif key == "vvel_cens":
            return self.vvel_cens

    def xymean(self, var):
        """mean over x and y"""

        return np.mean(var, axis=(1, 2))  # d ims [time, z]

    def ytmean(self, var):
        """mean over time and y"""

        return np.mean(var, axis=(0, 1))  # dims [x, z]

    def xytmean(self, var):
        """mean over time, x and y"""

        return np.mean(var, axis=(0, 1, 2))  # dims [z]


def thermovar_from_binary(var, thermofiles, shape, ntime, ndims, dtype, isprint=True):
    filename = f"{thermofiles.stem}_{var}{thermofiles.suffix}"
    filename = thermofiles.parent / filename
    data, ndata = readbinary(filename, isprint=isprint)

    if ndata != int(np.prod(shape)):
        err = (
            str(ndata)
            + " is incorrect data length for "
            + var
            + " defined for "
            + str(ntime)
            + " timesteps"
            + " on grid with dims = "
            + str(ndims)
        )
        raise ValueError(err)
    else:
        data = np.reshape(np.asarray(data, dtype=dtype), shape)

    return data


def read_dimless_thermodynamics_binary(thermofiles, ndims, ntime, nspacedims):
    # expected lengths of data defined on gridbox centres or faces
    cen = [ntime, int(np.prod(ndims))]
    zface = [ntime, int(ndims[2] * ndims[1] * (ndims[0] + 1))]
    xface = [ntime, int(ndims[2] * (ndims[1] + 1) * ndims[0])]
    yface = [ntime, int((ndims[2] + 1) * ndims[1] * ndims[0])]

    thermodata = {}

    vars = ["press", "temp", "qvap", "qcond"]
    datatypes = [np.double] * 4
    for v, var in enumerate(vars):
        thermodata[var] = thermovar_from_binary(
            var, thermofiles, cen, ntime, ndims, datatypes[v], isprint=False
        )

    datatypes = [np.double] * 3
    if nspacedims >= 1:
        thermodata["wvel"] = thermovar_from_binary(
            "wvel", thermofiles, zface, ntime, ndims, datatypes[0], isprint=False
        )
        if nspacedims >= 2:
            thermodata["uvel"] = thermovar_from_binary(
                "uvel", thermofiles, xface, ntime, ndims, datatypes[1], isprint=False
            )
            if nspacedims >= 3:
                thermodata["vvel"] = thermovar_from_binary(
                    "vvel",
                    thermofiles,
                    yface,
                    ntime,
                    ndims,
                    datatypes[2],
                    isprint=False,
                )

    return thermodata


def get_thermodynamics_from_thermofiles(
    thermofiles, ndims, inputs=False, constants_filename="", config_filename=""
):
    if not inputs:
        inputs = thermoinputsdict(config_filename, constants_filename)

    thermodata = read_dimless_thermodynamics_binary(
        thermofiles, ndims, inputs["ntime"], inputs["nspacedims"]
    )  # dimensionless data [time, gridboxes]
    thermodata = ThermoOnGrid(thermodata, inputs, ndims)

    return thermodata  # data with units in 4D arrays with dims [time, y, x, z]


def plot_thermodynamics(
    constants_filename,
    config_filename,
    grid_filename,
    thermofiles,
    savefig=False,
    savefigpath=None,
    savelabel="",
):
    plt.rcParams.update({"font.size": 14})

    inputs = thermoinputsdict(config_filename, constants_filename)
    gbxbounds, ndims = rgrid.read_dimless_gbxboundaries_binary(
        grid_filename, COORD0=inputs["COORD0"], return_ndims=True, isprint=False
    )
    xyzhalf = rgrid.halfcoords_from_gbxbounds(gbxbounds, isprint=False)  # [m]
    zhalf, xhalf, yhalf = [half / 1000 for half in xyzhalf]  # convery [m] to [km]
    zfull, xfull, yfull = rgrid.fullcell_fromhalfcoords(zhalf, xhalf, yhalf)  # [m]

    thermodata = get_thermodynamics_from_thermofiles(thermofiles, ndims, inputs=inputs)

    plot_1dprofiles(
        zfull,
        thermodata,
        inputs["Mr_ratio"],
        inputs["RGAS_DRY"],
        inputs["CP_DRY"],
        savefigpath,
        savefig,
        savelabel=savelabel,
    )

    if inputs["nspacedims"] > 1:
        xxh, zzh = np.meshgrid(xhalf, zhalf, indexing="ij")  # dims [xdims, zdims]
        xxf, zzf = np.meshgrid(xfull, zfull, indexing="ij")  # dims [xdims, zdims]
        plot_2dcolormaps(
            zzh,
            xxh,
            zzf,
            xxf,
            thermodata,
            inputs,
            savefigpath,
            savefig,
            savelabel=savelabel,
        )
        plot_2dwindfield(
            zzh,
            xxh,
            zzf,
            xxf,
            thermodata["wvel_cens"],
            thermodata["uvel_cens"],
            savefigpath,
            savefig,
            savelabel=savelabel,
        )


def try1dplot(ax, nplots, data, zfull, label):
    try:
        ax.plot(data, zfull, marker="x")  # (fails for 0D model)
    except ValueError as e:
        print("Plotting failed with ValueError:", e, " using ax.scatter instead.")
        ax.scatter(data, zfull, marker="x")

    ax.set_xlabel(label)

    return nplots + 1


def plot_1dthermodynamics(axs, n, zfull, thermodata, Mr_ratio, RGAS_DRY, CP_DRY):
    vars = ["press", "temp", "qvap", "qcond"]
    units = [" /Pa", " /K", "", ""]

    for var, unit in zip(vars, units):
        if var in thermodata.vars:
            profalltime = thermodata.xymean(thermodata[var])  # 1d profile at all times
            label = var + unit
            n = try1dplot(axs[n], n, profalltime.T, zfull[None, :].T, label)

    pressxy = thermodata.xymean(thermodata.press)
    tempxy = thermodata.xymean(thermodata.temp)
    qvapxy = thermodata.xymean(thermodata.qvap)
    supersat = relative_humidity(pressxy, tempxy, qvapxy, Mr_ratio)[1]
    label = "supersaturation"
    n = try1dplot(axs[n], n, supersat.T, zfull[None, :].T, label)

    press_ref = pressxy[0, 0]
    theta = potential_temperature(pressxy, tempxy, press_ref, RGAS_DRY, CP_DRY)
    label = "\u03F4 /K"
    n = try1dplot(axs[n], n, theta.T, zfull[None, :].T, label)

    return n


def plot_1dwindprofiles(axs, n, zfull, thermodata):
    vars = ["wvel_cens", "uvel_cens", "vvel_cens"]
    units = [" /ms$^{-1}$"] * 3

    for var, unit in zip(vars, units):
        if var in thermodata.vars:
            profalltime = thermodata.xymean(thermodata[var])  # 1d profile at all times
            label = var + unit
            n = try1dplot(axs[n], n, profalltime.T, zfull[None, :].T, label)

    return n


def plot_1dprofiles(
    zfull, thermodata, Mr_ratio, RGAS_DRY, CP_DRY, savefigpath, savefig, savelabel=""
):
    fig, axs = plt.subplots(nrows=3, ncols=3, figsize=(16, 8))
    axs = axs.flatten()

    nplots = plot_1dthermodynamics(
        axs, 0, zfull, thermodata, Mr_ratio, RGAS_DRY, CP_DRY
    )
    nplots = plot_1dwindprofiles(axs, nplots, zfull, thermodata)

    for a in range(nplots, len(axs), 1):
        axs[a].remove()  # delete unused axes
    for ax in [axs[0], axs[3], axs[6]]:
        ax.set_ylabel("z /km")

    fig.tight_layout()

    if savefig:
        savename = savefigpath / Path(f"thermo1dalltimeprofiles{savelabel}.png")
        fig.savefig(
            savename,
            dpi=400,
            bbox_inches="tight",
            facecolor="w",
            format="png",
        )
        print("Figure .png saved as: " + str(savefigpath))


def plot_2dmean(ax, meanfunc, xxh, zzh, var, label, unit, norm, cmap):
    mean2d = meanfunc(var)  # avg over time and y axes
    if np.nanmin(mean2d) != np.nanmax(mean2d):
        pcm = ax.pcolormesh(xxh[:, :], zzh[:, :], mean2d, cmap=cmap, norm=norm)
        cb = plt.colorbar(pcm, ax=ax, location="top", label=label + unit)
        return mean2d, pcm, cb
    else:
        txt = label + " = {:.2f}".format(np.nanmin(mean2d)) + unit
        ax.text(0.5, 0.5, txt, ha="center")
        return np.nanmin(mean2d), False, False


def plot_2dcontour(ax, xxf, zzf, mean2d, contour, cb, cbticks):
    cb.ax.set_xticks(cbticks)
    cb.ax.set_xticklabels(["{:.3g}".format(c) for c in cbticks])

    try:
        ax.contour(
            xxf, zzf, mean2d, levels=[contour], linestyles=["--"], colors=["grey"]
        )
    except ValueError as e:
        print("Plotting failed with ValueError:", e, "WARNING: not plotting contour")

    cb.ax.plot([contour] * 2, [0, 1], color="grey", linewidth=0.95)


def relh_supersat_theta_colomaps(
    axs, zzh, xxh, zzf, xxf, thermodata, Mr_ratio, RGAS_DRY, CP_DRY
):
    relh, supersat = relative_humidity(
        thermodata.press, thermodata.temp, thermodata.qvap, Mr_ratio
    )
    relh = relh * 100  # convert relative humidity to %

    press_ref = thermodata.xymean(thermodata.press)[0, 0]
    theta = potential_temperature(
        thermodata.press, thermodata.temp, press_ref, RGAS_DRY, CP_DRY
    )

    vars = [relh, supersat, theta]
    labels = ["% relative humidity", "supersaturation", "\u03F4"]
    units = ["", "", " /K"]
    cmaps = ["Blues", "PuOr", "RdBu_r"]
    norms = [
        colors.CenteredNorm(vcenter=np.mean(relh)),
        colors.TwoSlopeNorm(vcenter=0.0),
        colors.CenteredNorm(vcenter=np.nanmin(theta)),
    ]
    contours = [1.0, 0.0, None]

    n = 0
    for var, label, unit, norm, cmap, cont in zip(
        vars, labels, units, norms, cmaps, contours
    ):
        mean2d, pcm, cb = plot_2dmean(
            axs[n], thermodata.ytmean, xxh, zzh, var, label, unit, norm, cmap
        )

        if pcm is not False:
            if n == 1:
                cbticks = [pcm.norm.vmin, pcm.norm.vcenter, pcm.norm.vmax]
            else:
                cbticks = np.linspace(pcm.norm.vmin, pcm.norm.vmax, 5)
            plot_2dcontour(axs[n], xxf, zzf, mean2d, cont, cb, cbticks)

        n += 1


def plot_2dcolormaps(
    zzh, xxh, zzf, xxf, thermodata, inputs, savefigpath, savefig, savelabel=""
):
    vars = ["press", "temp", "qvap", "qcond"]
    units = [" /Pa", " /K", "", ""]
    cmaps = ["PRGn", "RdBu_r", "BrBG", "BrBG"]
    cmapcens = [
        np.nanmin(thermodata.press),
        np.nanmin(thermodata.temp),
        np.mean(thermodata.qvap),
        0.0,
    ]

    fig, axs = plt.subplots(nrows=2, ncols=4, figsize=(13, 8))
    axs = axs.flatten()
    axs[-1].remove()

    n = 0
    for var, unit, cmapcen, cmap in zip(vars, units, cmapcens, cmaps):
        norm = colors.CenteredNorm(vcenter=cmapcen)
        plot_2dmean(
            axs[n], thermodata.ytmean, xxh, zzh, thermodata[var], var, unit, norm, cmap
        )
        n += 1

    relh_supersat_theta_colomaps(
        [axs[4], axs[5], axs[6]],
        zzh,
        xxh,
        zzf,
        xxf,
        thermodata,
        inputs["Mr_ratio"],
        inputs["RGAS_DRY"],
        inputs["CP_DRY"],
    )

    for ax in axs:
        ax.set_aspect("equal")

    axs[0].set_ylabel("z /km")
    axs[3].set_ylabel("z /km")
    for ax in axs[3:]:
        ax.set_xlabel("x /km")

    fig.tight_layout()

    if savefig:
        savename = savefigpath / Path(f"thermo2dmeanprofiles{savelabel}.png")
        fig.savefig(
            savename,
            dpi=400,
            bbox_inches="tight",
            facecolor="w",
            format="png",
        )
        print("Figure .png saved as: " + str(savename))


def plot_2dwindfield(
    zzh, xxh, zzf, xxf, wvel_cens, uvel_cens, savefigpath, savefig, savelabel=""
):
    wcen = np.mean(wvel_cens, axis=(0, 1))  # avg over y and time axes
    ucen = np.mean(uvel_cens, axis=(0, 1))
    norm = np.sqrt(wcen**2 + ucen**2)

    fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(9, 6))
    axs = axs.flatten()

    label = ["w", "u", "|wind velocity|"]
    units = " / m s$^{-1}$"
    cmaps = ["coolwarm"] * 3
    for v, vel in enumerate([wcen, ucen, norm]):
        pcm = axs[v].pcolormesh(xxh, zzh, vel, cmap=cmaps[v])
        plt.colorbar(pcm, ax=axs[v], location="top", label=label[v] + units)
    axs[2].quiver(xxf, zzf, ucen, wcen)

    axs[0].set_ylabel("z /km")
    for ax in axs:
        ax.set_xlabel("x /km")
        ax.set_aspect("equal")

    fig.tight_layout()

    if savefig:
        savename = savefigpath / Path(f"thermowindprofiles{savelabel}.png")
        fig.savefig(
            savename,
            dpi=400,
            bbox_inches="tight",
            facecolor="w",
            format="png",
        )
        print("Figure .png saved as: " + str(savename))
