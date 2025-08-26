"""
Copyright (c) 2024 MPI-M, Clara Bayley


----- CLEO -----
File: breakup.py
Project: boxmodelcollisions
Created Date: Friday 14th June 2024
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
Script generates input files, runs CLEO 0-D box model executables for
collisions with selected collision kernels with breakup (e.g. Low and Lists's)
to create data. Then plots results analgous to Shima et al. 2009 Fig. 2(b)
"""

# %%
### -------------------------------- IMPORTS ------------------------------- ###
import argparse
import shutil
import subprocess
import sys
from pathlib import Path
from ruamel.yaml import YAML

# %%
### --------------------------- PARSE ARGUMENTS ---------------------------- ###
parser = argparse.ArgumentParser()
parser.add_argument(
    "path2CLEO", type=Path, help="Absolute path to CLEO directory (for cleopy)"
)
parser.add_argument("path2build", type=Path, help="Absolute path to build directory")
parser.add_argument(
    "src_config_filename",
    type=Path,
    help="Absolute path to source configuration YAML file",
)
parser.add_argument(
    "--kernels",
    nargs="*",
    choices=["long", "lowlist", "szakallurbich", "testikstraub"],
    type=str,
    help="kernel examples to run",
)
parser.add_argument(
    "--do_inputfiles",
    action="store_true",  # default is False
    help="Generate initial condition binary files",
)
parser.add_argument(
    "--do_run_executable",
    action="store_true",  # default is False
    help="Run executable",
)
parser.add_argument(
    "--do_plot_results",
    action="store_true",  # default is False
    help="Plot results of example",
)
args = parser.parse_args()

# %%
### -------------------------- INPUT PARAMETERS ---------------------------- ###
### --- command line parsed arguments --- ###
path2CLEO = args.path2CLEO
path2build = args.path2build
src_config_filename = args.src_config_filename
kernels = args.kernels

### --- additional/derived arguments --- ###
tmppath = path2build / "tmp"
sharepath = path2build / "share"
binpath = path2build / "bin"
savefigpath = binpath

kernel_configs = {}  # kernel: [config_filename, config_params]

isfigures = [False, True]  # booleans for [showing, saving] initialisation figures

for k in kernels:
    cf = path2build / "tmp" / f"breakup_{k}_config.yaml"
    cp = {
        "constants_filename": str(path2CLEO / "libs" / "cleoconstants.hpp"),
        "grid_filename": str(sharepath / "breakup_dimlessGBxboundaries.dat"),
        "maxnsupers": 8192,
        "initsupers_filename": str(sharepath / f"breakup_{k}_dimlessSDsinit.dat"),
        "setup_filename": str(binpath / f"breakup_{k}_setup.txt"),
        "zarrbasedir": str(binpath / f"breakup_{k}_sol.zarr"),
    }
    kernel_configs[k] = [cf, cp]

executables = {
    "long": "longcolls",
    "lowlist": "lowlistcolls",
    "szakallurbich": "szakallurbichcolls",
    "testikstraub": "testikstraubcolls",
}


# %%
### ------------------------- FUNCTION DEFINITIONS ------------------------- ###
def inputfiles(
    path2CLEO,
    path2build,
    tmppath,
    sharepath,
    binpath,
    savefigpath,
    src_config_filename,
    config_filename,
    config_params,
    kernel,
    isfigures,
):
    from cleopy import editconfigfile

    ### --- ensure build, share and bin directories exist --- ###
    if path2CLEO == path2build:
        raise ValueError("build directory cannot be CLEO")
    path2build.mkdir(exist_ok=True)
    tmppath.mkdir(exist_ok=True)
    sharepath.mkdir(exist_ok=True)
    binpath.mkdir(exist_ok=True)
    if savefigpath is not None:
        savefigpath.mkdir(exist_ok=True)

    ### --- copy src_config_filename into tmp and edit parameters --- ###
    config_filename.unlink(missing_ok=True)  # delete any existing config
    shutil.copy(src_config_filename, config_filename)
    editconfigfile.edit_config_params(config_filename, config_params)

    ### --- delete any existing initial conditions --- ###
    yaml = YAML()
    with open(config_filename, "r") as file:
        config = yaml.load(file)
    Path(config["inputfiles"]["grid_filename"]).unlink(missing_ok=True)
    Path(config["initsupers"]["initsupers_filename"]).unlink(missing_ok=True)

    ### --- input binary files generation --- ###
    # equivalent to ``import breakup_inputfiles`` followed by
    # ``breakup_inputfiles.main(path2CLEO, path2build, ...)``
    inputfiles_script = (
        path2CLEO / "examples" / "boxmodelcollisions" / "breakup_inputfiles.py"
    )
    python = sys.executable
    cmd = [
        python,
        inputfiles_script,
        path2CLEO,
        path2build,
        config_filename,
        kernel,
    ]
    if isfigures[0]:
        cmd.append("--show_figures")
    if isfigures[1]:
        cmd.append("--save_figures")
        cmd.append(f"--savefigpath={savefigpath}")
    print(" ".join([str(c) for c in cmd]))
    subprocess.run(cmd, check=True)


def run_exectuable(executable, config_filename):
    ### --- delete any existing output dataset and setup files --- ###
    yaml = YAML()
    with open(config_filename, "r") as file:
        config = yaml.load(file)
    Path(config["outputdata"]["setup_filename"]).unlink(missing_ok=True)
    shutil.rmtree(Path(config["outputdata"]["zarrbasedir"]), ignore_errors=True)

    ### --- run exectuable with given config file --- ###
    cmd = [executable, config_filename]
    print(" ".join([str(c) for c in cmd]))
    subprocess.run(cmd, check=True)


def plot_results(path2CLEO, savefigpath, kernels, config_filenames):
    plotting_script = (
        path2CLEO / "examples" / "boxmodelcollisions" / "breakup_plotting.py"
    )
    python = sys.executable

    setupfiles, datasets = [], []
    for config_filename in config_filenames:
        yaml = YAML()
        with open(config_filename, "r") as file:
            config = yaml.load(file)
        setupfiles.append(Path(config["outputdata"]["setup_filename"]))
        datasets.append(Path(config["outputdata"]["zarrbasedir"]))
    grid_filename = Path(
        config["inputfiles"]["grid_filename"]
    )  # same grid for all datasets

    # equivalent to ``import breakup_plotting`` followed by
    # ``breakup_plotting.main(path2CLEO, savefigpath, ...)``
    cmd = [
        python,
        plotting_script,
        f"--path2CLEO={path2CLEO}",
        f"--savefigpath={savefigpath}",
        f"--grid_filename={grid_filename}",
    ]
    cmd += ["--kernels"] + kernels
    cmd += ["--setupfiles"] + setupfiles
    cmd += ["--datasets"] + datasets
    print(" ".join([str(c) for c in cmd]))
    subprocess.run(cmd, check=True)


# %%
### --------------------- RUN EXAMPLE FOR EACH KERNEL ---------------------- ###
for kernel, [config_filename, config_params] in kernel_configs.items():
    if args.do_inputfiles:
        inputfiles(
            path2CLEO,
            path2build,
            tmppath,
            sharepath,
            binpath,
            savefigpath,
            src_config_filename,
            config_filename,
            config_params,
            kernel,
            isfigures,
        )

    if args.do_run_executable:
        executable = (
            path2build / "examples" / "boxmodelcollisions" / "src" / executables[kernel]
        )
        run_exectuable(executable, config_filename)

if args.do_plot_results:
    kernels2plot = list(kernel_configs.keys())
    config_filenames2plot = list([cfgs[0] for cfgs in kernel_configs.values()])
    plot_results(path2CLEO, savefigpath, kernels2plot, config_filenames2plot)
