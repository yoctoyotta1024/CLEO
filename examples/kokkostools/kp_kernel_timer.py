"""
Copyright (c) 2025 MPI-M, Clara Bayley


----- CLEO -----
File: kp_kernel_timer.py
Project: kokkostools
Created Date: Monday 25th August 2025
Author: Clara Bayley (CB)
Additional Contributors:
-----
License: BSD 3-Clause "New" or "Revised" License
https://opensource.org/licenses/BSD-3-Clause
-----
File Description:
Class for using the Kokkos Tools Kernel Timer within Python
"""


# %%
### -------------------------------- IMPORTS ------------------------------- ###
import glob
import os
import subprocess
from pathlib import Path


# %%
### --------------------- KOKKOS PROFILER DEFINITION ----------------------- ###
class KpKernelTimer:
    def __init__(self, path2kokkostools: Path):
        print("---------------------------------------------------")
        print(f"Using Kokkos Tools in: {path2kokkostools}")
        assert path2kokkostools.is_dir(), "path2kokkostools doesn't exist"

        kokkostools_lib = glob.glob(os.path.join(path2kokkostools, "lib*/"))[0]
        self.kokkostools_lib = Path(kokkostools_lib)

        kp_kernel_timer = glob.glob(
            os.path.join(self.kokkostools_lib, "libkp_kernel_timer.*")
        )[0]
        os.environ["KOKKOS_TOOLS_LIBS"] = str(kp_kernel_timer)

        self.kp_reader = path2kokkostools / "bin" / "kp_reader"

        print("Using Kokkos Tool libs:", self.kokkostools_lib)
        print("Using Kokkos Profiling Tool:", os.environ["KOKKOS_TOOLS_LIBS"])
        print("Using Kokkos Tool Reader:", self.kp_reader)
        print("---------------------------------------------------")

        assert Path(
            os.environ["KOKKOS_TOOLS_LIBS"]
        ).exists(), "kp_kernel_timer lib doesn't exist"
        assert self.kp_reader.exists(), "kp_reader doesn't exist"

    def postprocess(self, data_filespath: Path, txt_filename: Path):
        assert (
            txt_filename.parent
        ).is_dir(), "profiler data output directory doesn't exist"
        assert txt_filename.suffix == ".txt", "profiler data output must be .txt file"

        # Add kokkostools_lib to LD_LIBRARY_PATH
        assert self.kokkostools_lib.is_dir(), "kokkostools_lib doesn't exist"
        ld_lib_path = os.environ.get("LD_LIBRARY_PATH", "")  # linux systems
        os.environ["LD_LIBRARY_PATH"] = f"{self.kokkostools_lib}:{ld_lib_path}"
        dyld_lib_path = os.environ.get("DYLD_LIBRARY_PATH", "")  # macOS systems
        os.environ["DYLD_LIBRARY_PATH"] = f"{self.kokkostools_lib}:{dyld_lib_path}"

        # Use glob to find all .dat data files in the specified directory
        datfiles = glob.glob(os.path.join(data_filespath, "*.dat"))

        # Post-process each .dat file and write to txt_filename_{f}.txt
        for f, datafile in enumerate(datfiles):
            datafile = Path(datafile)
            txt_filename_run = (
                f"{txt_filename.stem}_{f}_{datafile.stem}{txt_filename.suffix}"
            )
            txt_filename_run = txt_filename.parent / txt_filename_run
            cmd = [str(self.kp_reader), str(datafile)]
            with open(txt_filename_run, "w") as wfile:
                subprocess.run(cmd, stdout=wfile, stderr=subprocess.STDOUT, check=True)
