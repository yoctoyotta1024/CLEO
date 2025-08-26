"""
Copyright (c) 2024 MPI-M, Clara Bayley


----- CLEO -----
File: shima2009.py
Project: boxmodelcollisions
Created Date: Friday 17th November 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
Script generates input files, runs CLEO 0-D box model executables for collisions
with selected collision kernels (e.g. Golovin's or Long's) to create data.
Then plots results comparable to Shima et al. 2009 Fig. 2
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
    choices=["golovin", "long1", "long2"],
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

if "golovin" in kernels:
    k = "golovin"
    cf = path2build / "tmp" / f"shima2009_{k}_config.yaml"
    cp = {
        "constants_filename": str(path2CLEO / "libs" / "cleoconstants.hpp"),
        "grid_filename": str(sharepath / "shima2009_dimlessGBxboundaries.dat"),
        "COLLTSTEP": 1,
        "maxnsupers": 4096,
        "initsupers_filename": str(sharepath / f"shima2009_{k}_dimlessSDsinit.dat"),
        "setup_filename": str(binpath / f"shima2009_{k}_setup.txt"),
        "zarrbasedir": str(binpath / f"shima2009_{k}_sol.zarr"),
    }
    kernel_configs[k] = [cf, cp]

if "long1" in kernels:
    k = "long1"
    cf = path2build / "tmp" / f"shima2009_{k}_config.yaml"
    cp = {
        "constants_filename": str(path2CLEO / "libs" / "cleoconstants.hpp"),
        "grid_filename": str(sharepath / "shima2009_dimlessGBxboundaries.dat"),
        "COLLTSTEP": 1,
        "maxnsupers": 4096,
        "initsupers_filename": str(sharepath / f"shima2009_{k}_dimlessSDsinit.dat"),
        "setup_filename": str(binpath / f"shima2009_{k}_setup.txt"),
        "zarrbasedir": str(binpath / f"shima2009_{k}_sol.zarr"),
    }
    kernel_configs[k] = [cf, cp]

if "long2" in kernels:
    k = "long2"
    cf = path2build / "tmp" / f"shima2009_{k}_config.yaml"
    cp = {
        "constants_filename": str(path2CLEO / "libs" / "cleoconstants.hpp"),
        "grid_filename": str(sharepath / "shima2009_dimlessGBxboundaries.dat"),
        "COLLTSTEP": 0.1,
        "maxnsupers": 8192,
        "initsupers_filename": str(sharepath / f"shima2009_{k}_dimlessSDsinit.dat"),
        "setup_filename": str(binpath / f"shima2009_{k}_setup.txt"),
        "zarrbasedir": str(binpath / f"shima2009_{k}_sol.zarr"),
    }
    kernel_configs[k] = [cf, cp]

executables = {
    "golovin": "golcolls",
    "long1": "longcolls",
    "long2": "longcolls",
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
    # equivalent to ``import shima2009_inputfiles`` followed by
    # ``shima2009_inputfiles.main(path2CLEO, path2build, ...)``
    inputfiles_script = (
        path2CLEO / "examples" / "boxmodelcollisions" / "shima2009_inputfiles.py"
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


def plot_results(path2CLEO, config_filename, savefigpath, kernel):
    plotting_script = (
        path2CLEO / "examples" / "boxmodelcollisions" / "shima2009_plotting.py"
    )
    python = sys.executable

    yaml = YAML()
    with open(config_filename, "r") as file:
        config = yaml.load(file)
    grid_filename = Path(config["inputfiles"]["grid_filename"])
    setupfile = Path(config["outputdata"]["setup_filename"])
    dataset = Path(config["outputdata"]["zarrbasedir"])

    # equivalent to ``import shima2009_plotting`` followed by
    # ``shima2009_plotting.main(path2CLEO, savefigpath, ...)``
    cmd = [
        python,
        plotting_script,
        f"--path2CLEO={path2CLEO}",
        f"--savefigpath={savefigpath}",
        f"--grid_filename={grid_filename}",
        f"--setupfile={setupfile}",
        f"--dataset={dataset}",
        f"--kernel={kernel}",
    ]
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
        plot_results(path2CLEO, config_filename, savefigpath, kernel)
