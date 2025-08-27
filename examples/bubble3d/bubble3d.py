"""
Copyright (c) 2024 MPI-M, Clara Bayley


----- CLEO -----
File: bubble3d.py
Project: bubble3d
Created Date: Friday 17th November 2023
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
Script generates input files, then runs CLEO executable "bubble3d" to
piggyback ICON bubble test case
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
    help="Run example executable(s)",
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

config_filename = path2build / "tmp" / "bubble3d_config.yaml"
config_params = {
    "constants_filename": str(path2CLEO / "libs" / "cleoconstants.hpp"),
    "grid_filename": str(sharepath / "bubble3d_dimlessGBxboundaries.dat"),
    "initsupers_filename": str(sharepath / "bubble3d_dimlessSDsinit.dat"),
    "setup_filename": str(binpath / "bubble3d_setup.txt"),
    "zarrbasedir": str(binpath / "bubble3d_sol.zarr"),
    "orginal_icon_grid_file": "/work/bm1183/m300950/icon/build/experiments/aes_bubble/aes_bubble_atm_cgrid_ml.nc",
    "orginal_icon_data_file": "/work/bm1183/m300950/icon/build/experiments/aes_bubble/aes_bubble_atm_3d_ml_20080801T000000Z.nc",
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
    isfigures,
):
    from cleopy import editconfigfile

    ### --- ensure build, tmp, share and bin and savefigpath directories exist --- ###
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
    # equivalent to ``import bubble3d_inputfiles`` followed by
    # ``bubble3d_inputfiles.main(path2CLEO, path2build, ...)``
    inputfiles_script = path2CLEO / "examples" / "bubble3d" / "bubble3d_inputfiles.py"
    python = sys.executable
    cmd = [python, inputfiles_script, path2CLEO, path2build, config_filename]
    if isfigures[0]:
        cmd.append("--show_figures")
    if isfigures[1]:
        cmd.append("--save_figures")
        cmd.append(f"--savefigpath={savefigpath}")
    print(" ".join([str(c) for c in cmd]))
    subprocess.run(cmd, check=True)


def run_exectuable(path2CLEO, path2build, config_filename):
    ### --- delete any existing output dataset and setup files --- ###
    yaml = YAML()
    with open(config_filename, "r") as file:
        config = yaml.load(file)
    Path(config["outputdata"]["setup_filename"]).unlink(missing_ok=True)
    shutil.rmtree(Path(config["outputdata"]["zarrbasedir"]), ignore_errors=True)

    ### --- run exectuable with given config file --- ###
    cleoproc = str(path2build / "examples" / "bubble3d" / "src" / "bubble3d")
    cleoproc_args = str(config_filename)
    print("CLEO Executable: " + cleoproc)
    print("CLEO Config file: " + cleoproc_args)

    python = sys.executable
    pythonproc = str(path2CLEO / "examples" / "bubble3d" / "yac_bubble_data_reader.py")
    pythonproc_args = [str(path2build), str(config_filename)]
    print("YAC script: " + pythonproc)
    print("YAC arguments: " + " ".join(pythonproc_args))

    cmd = [
        "mpiexec",
        "-n",
        "1",
        cleoproc,
        cleoproc_args,
        ":",
        "-n",
        "1",
        python,
        pythonproc,
    ] + pythonproc_args
    print(" ".join([str(c) for c in cmd]))
    subprocess.run(cmd, check=True)


def plot_results(path2CLEO, config_filename, savefigpath):
    plotting_script = path2CLEO / "examples" / "bubble3d" / "bubble3d_plotting.py"
    python = sys.executable

    yaml = YAML()
    with open(config_filename, "r") as file:
        config = yaml.load(file)
    grid_filename = Path(config["inputfiles"]["grid_filename"])
    setupfile = Path(config["outputdata"]["setup_filename"])
    dataset = Path(config["outputdata"]["zarrbasedir"])

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
        isfigures,
    )

if args.do_run_executable:
    run_exectuable(path2CLEO, path2build, config_filename)

if args.do_plot_results:
    plot_results(path2CLEO, config_filename, savefigpath)
