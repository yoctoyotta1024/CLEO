import glob
import re
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

# is_finddata = True
# is_plot = False
is_finddata = False
is_plot = True

if is_finddata:
    ntasks = 4
    time = []
    for file_path in glob.glob(
        f"/home/m/m300950/CLEO/build_fromfile/bin/ntasks{ntasks}/mpiscaling_out.*.out"
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
        1: [564.5],
        4: [594.74, 594.62, 594.67, 594.67],
        8: [563.09, 562.95, 562.95, 562.95, 562.95, 562.92, 562.92, 562.93],
        16: [
            543.04,
            543.04,
            543.04,
            543.04,
            543.04,
            543.04,
            543.04,
            543.04,
            543.02,
            543.02,
            543.02,
            543.02,
            543.02,
            543.02,
            543.02,
            543.17,
        ],
        32: [
            540.66,
            540.66,
            540.66,
            540.65,
            540.66,
            540.65,
            540.65,
            540.65,
            540.65,
            540.65,
            540.65,
            540.66,
            540.66,
            540.65,
            540.65,
            540.66,
            540.66,
            540.66,
            540.66,
            540.65,
            540.66,
            540.66,
            540.65,
            540.66,
            540.65,
            540.65,
            540.65,
            540.66,
            540.65,
            540.66,
            540.66,
            540.78,
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
