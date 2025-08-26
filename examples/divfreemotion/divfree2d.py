"""
Copyright (c) 2024 MPI-M, Clara Bayley


----- CLEO -----
File: divfree2d.py
Project: divfreemotion
Created Date: Friday 17th November 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
Script generates input files, then runs CLEO executable "divfree2d" to create the
data for example of divergence free motion of superdroplets in a 2-D divergence
free flow field.
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

### --- additional/derived arguments --- ###
tmppath = path2build / "tmp"
sharepath = path2build / "share"
binpath = path2build / "bin"
savefigpath = binpath

config_filename = path2build / "tmp" / "divfree2d_config.yaml"
thermofiles = sharepath / "divfree2d_dimlessthermo.dat"
config_params = {
    "constants_filename": str(path2CLEO / "libs" / "cleoconstants.hpp"),
    "grid_filename": str(sharepath / "divfree2d_dimlessGBxboundaries.dat"),
    "initsupers_filename": str(sharepath / "divfree2d_dimlessSDsinit.dat"),
    "setup_filename": str(binpath / "divfree2d_setup.txt"),
    "zarrbasedir": str(binpath / "divfree2d_sol.zarr"),
}

isfigures = [False, True]  # booleans for [showing, saving] initialisation figures


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
    thermofiles,
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

    ### --- add names of thermofiles to config_params --- ###
    for var in ["press", "temp", "qvap", "qcond", "wvel", "uvel"]:
        config_params[var] = str(
            thermofiles.parent / Path(f"{thermofiles.stem}_{var}{thermofiles.suffix}")
        )

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
    all_thermofiles = thermofiles.parent.glob(
        f"{thermofiles.stem}*{thermofiles.suffix}"
    )
    for file in all_thermofiles:
        file.unlink(missing_ok=True)

    ### --- input binary files generation --- ###
    # equivalent to ``import divfree2d_inputfiles`` followed by
    # ``divfree2d_inputfiles.main(path2CLEO, path2build, ...)``
    inputfiles_script = (
        path2CLEO / "examples" / "divfreemotion" / "divfree2d_inputfiles.py"
    )
    python = sys.executable
    cmd = [
        python,
        inputfiles_script,
        path2CLEO,
        path2build,
        config_filename,
        thermofiles,
    ]
    if isfigures[0]:
        cmd.append("--show_figures")
    if isfigures[1]:
        cmd.append("--save_figures")
        cmd.append(f"--savefigpath={savefigpath}")
    print(" ".join([str(c) for c in cmd]))
    subprocess.run(cmd, check=True)


def run_exectuable(path2build, config_filename):
    ### --- delete any existing output dataset and setup files --- ###
    yaml = YAML()
    with open(config_filename, "r") as file:
        config = yaml.load(file)
    Path(config["outputdata"]["setup_filename"]).unlink(missing_ok=True)
    shutil.rmtree(Path(config["outputdata"]["zarrbasedir"]), ignore_errors=True)

    ### --- run exectuable with given config file --- ###
    executable = path2build / "examples" / "divfreemotion" / "src" / "divfree2d"
    cmd = [executable, config_filename]
    print(" ".join([str(c) for c in cmd]))
    subprocess.run(cmd, check=True)


def plot_results(path2CLEO, config_filename, savefigpath):
    plotting_script = path2CLEO / "examples" / "divfreemotion" / "divfree2d_plotting.py"
    python = sys.executable

    yaml = YAML()
    with open(config_filename, "r") as file:
        config = yaml.load(file)
    grid_filename = Path(config["inputfiles"]["grid_filename"])
    setupfile = Path(config["outputdata"]["setup_filename"])
    dataset = Path(config["outputdata"]["zarrbasedir"])

    # equivalent to ``import divfree2d_plotting`` followed by
    # ``divfree2d_plotting.main(path2CLEO, savefigpath, ...)``
    cmd = [
        python,
        plotting_script,
        f"--path2CLEO={path2CLEO}",
        f"--savefigpath={savefigpath}",
        f"--grid_filename={grid_filename}",
        f"--setupfile={setupfile}",
        f"--dataset={dataset}",
    ]
    print(" ".join([str(c) for c in cmd]))
    subprocess.run(cmd, check=True)


# %%
### ----------------------------- RUN EXAMPLE ------------------------------ ###
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
        thermofiles,
        isfigures,
    )

if args.do_run_executable:
    run_exectuable(path2build, config_filename)

if args.do_plot_results:
    plot_results(path2CLEO, config_filename, savefigpath)
