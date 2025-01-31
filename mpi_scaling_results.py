import glob
import re
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

is_finddata = False
is_plot = True

if is_finddata:
    time = []
    for file_path in glob.glob(
        "/home/m/m300950/CLEO/build_fromfile/bin/mpiscaling_out.*.out"
    ):
        with open(file_path, "r") as file:
            lines = file.readlines()

            for line in lines:
                if "Total Program Duration" in line:
                    match = re.search(r"(\d+\.\d+e\+\d+)", line.strip())
                    time.append(float(match.group(1)))
    print("time: ", time)


if is_plot:
    raw_data = {
        1: [1.0852e02, 1.0844e02],
        4: [1.9338e02, 1.9349e02, 1.9338e02, 1.9337e02],
        8: [
            1.7531e02,
            1.7531e02,
            1.7531e02,
            1.7531e02,
            1.7531e02,
            1.7530e02,
            1.7530e02,
            1.7544e02,
        ],
        16: [1.7611e02] * 12 + [1.7612e02] * 3 + [1.7624e02],
        32: [
            176.55,
            176.55,
            176.55,
            176.55,
            176.55,
            176.55,
            176.55,
            176.55,
            176.54,
            176.54,
            176.55,
            176.54,
            176.55,
            176.55,
            176.55,
            176.54,
            176.54,
            176.55,
            176.55,
            176.55,
            176.55,
            176.55,
            176.54,
            176.54,
            176.54,
            176.55,
            176.55,
            176.55,
            176.54,
            176.55,
            176.54,
            176.68,
        ],
    }

    ntasks = list(raw_data.keys())
    time = xr.DataArray(
        [np.mean(raw_data[v]) for v in raw_data.keys()],
        dims=["ntasks"],
        attrs={"units": "s", "long-name": "wall-clock time"},
    )

    ds = xr.Dataset(data_vars={"time": time}, coords={"ntasks": (ntasks)})
    print(ds)

    fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(10, 6), sharex=True)
    fig.suptitle("Strong Scaling of Test Case for CLEO's MPI Domain Decomposition")
    axs[1].set_xlabel("number of MPI processes")

    ds["time"].plot(ax=axs[0], marker="o", linestyle="-")
    axs[0].set_ylabel(f"total wall-clock time / {ds.time.units}")

    speedup = ds["time"].sel(ntasks=1) / ds["time"]
    speedup.plot(ax=axs[1], marker="o", linestyle="-")
    axs[1].set_ylabel("speedup")

    efficiency = ds["time"].sel(ntasks=1) / ds["time"] / ds["ntasks"]
    efficiency.plot(ax=axs[2], marker="o", linestyle="-")
    axs[2].set_ylabel("efficiency")

    for ax in axs:
        ax.spines[["right", "top"]].set_visible(False)

    plt.tight_layout()
    plt.savefig("./strong_scaling_wallclock_time.png", bbox_inches="tight")
