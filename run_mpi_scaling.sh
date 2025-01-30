#!/bin/bash
#SBATCH --job-name=mpiscaling
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=16
#SBATCH --mem=30G
#SBATCH --time=00:10:00
#SBATCH --mail-user=clara.bayley@mpimet.mpg.de
#SBATCH --mail-type=FAIL
#SBATCH --account=bm1183
#SBATCH --output=./mpiscaling_out.%j.out
#SBATCH --error=./mpiscaling_err.%j.out

### ---------------------------------------------------- ###
### Useful Info for Profiling An Executable Using Kokkos ###
### ---------------------------------------------------- ###
# E.g. for Kokkos tools installed in /work/bm1183/m300950/kokkos_tools_lib/:
# A) see tool libraries installed in /work/bm1183/m300950/kokkos_tools_lib/lib64/
# B) export required tool library, e.g.
#     e.g. export KOKKOS_TOOLS_LIBS=/work/bm1183/m300950/kokkos_tools_lib/lib64/libkp_kernel_timer.so
#      or  export KOKKOS_TOOLS_LIBS=/work/bm1183/m300950/kokkos_tools_lib/lib64/libkp_space_time_stack.so
# C) run executable ./[exec].exe (kokkos initialise loads dynamic library pointers)
# D) read *.dat output
#     e.g. with kp reader
#          export LD_LIBRARY_PATH=/work/bm1183/m300950/kokkos_tools_lib/lib64/:$LD_LIBRARY_PATH
#          /work/bm1183/m300950/kokkos_tools_lib/bin/kp_reader *.dat > ./bin/kp_kernel_timer.txt
#     or pipe kp_space_time_stack output durign runtime:
#          ./[exec].exe > runtime_output.txt
#
# + see useful debugging tool to find where program crashed (e.g. inside kernel):
#     export KOKKOS_TOOLS_LIBS=/work/bm1183/m300950/kokkos_tools_lib/lib64/libkp_kernel_logger.so
### ---------------------------------------------------- ###

### ---------------------------------------------------- ###
### ------------------ Input Parameters ---------------- ###
### ------ You MUST edit these lines to set your ------- ###
### ---- build type, directories, the executable(s) ---- ###
### -------- to compile, and your python script -------- ###
### ---------------------------------------------------- ###
actions=("$@")

ntasks=8 ### Useful: https://slurm-binding.dkrz.de/singlenodebinding.html

buildtype="openmp"
path2CLEO=${HOME}/CLEO/
path2build=${HOME}/CLEO/build_fromfile/
enableyac=false
executables="fromfile"

pythonscript=${path2CLEO}/examples/fromfile/fromfile.py
configfile_og=${path2CLEO}/examples/fromfile/src/config/fromfile_config.yaml
configfile=${path2build}/tmp/fromfile_ntasks${ntasks}_config.yaml
cleoenv=/work/bm1183/m300950/bin/envs/cleoenv
python=${cleoenv}/bin/python3
enabledebug=false
make_clean=false
yacyaxtroot=/work/bm1183/m300950/yacyaxt
stacksize_limit=204800 # ulimit -s [stacksize_limit] (kB)
compilername=intel
### ---------------------------------------------------- ###
### ---------------------------------------------------- ###
### ---------------------------------------------------- ###

### -------------------- print inputs ------------------ ###
if [ ${#actions[@]} -eq 0 ]; then
    echo "Please specify action(s) out of: compile, inputfiles, run, plot, compare and delete"
    exit 1
fi
echo "----- Running Example -----"
echo "buildtype = ${buildtype}"
echo "compilername = ${compilername}"
echo "path2CLEO = ${path2CLEO}"
echo "path2build = ${path2build}"
echo "enableyac = ${enableyac}"
echo "executables = ${executables}"
echo "pythonscript = ${pythonscript}"
echo "---------------------------"
### ---------------------------------------------------- ###

### --------------- build and compile CLEO ------------- ###
if [[ " ${actions[@]} " =~ " compile " ]]; then
  cmd="${path2CLEO}/scripts/build_compile_cleo.sh \
    ${buildtype}
    ${compilername}
    ${path2CLEO}
    ${path2build}
    "\"${executables}\""
    ${enabledebug}
    ${enableyac}
    ${yacyaxtroot}
    ${make_clean}"
  echo ${cmd}
  eval ${cmd}
fi
### ---------------------------------------------------- ###

### ------ Gen Input Files through Python script ------- ###
if [[ " ${actions[@]} " =~ " inputfiles " ]]; then
  cp ${configfile_og} ${configfile}
  ${python} ${pythonscript} ${path2CLEO} ${path2build} ${configfile} \
    --do_inputfiles=TRUE --do_run_executable=FALSE --do_plot_results=FALSE --ntasks=${ntasks}
fi
### ---------------------------------------------------- ###

### -------- Run Model with Kokkos Kernel Timer -------- ###
if [[ " ${actions[@]} " =~ " run " ]]; then
  cd ${path2build} && pwd
  rm -rf ${path2build}/bin/ntasks${ntasks}/fromfile_sol.zarr

  export CLEO_PATH2CLEO=${path2CLEO}
  export CLEO_BUILDTYPE=${buildtype}
  export CLEO_ENABLEYAC=${enableyac}
  export KOKKOS_TOOLS_LIBS=/work/bm1183/m300950/kokkos_tools_lib/lib64/libkp_kernel_timer.so
  export LD_LIBRARY_PATH=/work/bm1183/m300950/kokkos_tools_lib/lib64/:$LD_LIBRARY_PATH
  source ${path2CLEO}/scripts/bash/src/runtime_settings.sh ${stacksize_limit}

  exec="srun --ntasks=${ntasks} ${path2build}/examples/fromfile/src/fromfile ${configfile}"
  echo ${exec}
  eval ${exec}

  /work/bm1183/m300950/kokkos_tools_lib/bin/kp_reader *.dat > ./bin/ntasks${ntasks}/kp_kerneltimer_ntasks${ntasks}.txt
fi
### ---------------------------------------------------- ###

### ------- Plot Results through Python script --------- ###
if [[ " ${actions[@]} " =~ " plot " ]]; then
  ${python} ${pythonscript} ${path2CLEO} ${path2build} ${configfile} \
    --do_inputfiles=FALSE --do_run_executable=FALSE --do_plot_results=TRUE --ntasks=${ntasks}
fi
### ---------------------------------------------------- ###

### --------- Check Results (after handling!) ---------- ###
if [[ " ${actions[@]} " =~ " compare " ]]; then
  cd ${path2build} && pwd
  if [ -s "./diffs" ]; then
    rm diffs
  fi
  for file in $(ls bin/ntasks8/fromfile_sol.zarr/); do
    for filename in $(ls bin/ntasks8/fromfile_sol.zarr/$file); do
      hexdump bin/ntasks8/fromfile_sol.zarr/$file/$filename > original
      hexdump bin/ntasks${ntasks}/fromfile_sol.zarr/$file/$filename > new
      diff original new >> diffs
    done
  done;
  rm original new

  if [ -s "./diffs" ]; then
    echo "ERROR: Run with ${ntasks} processes has different results than control run"
    exit 1;
  fi
fi
### ---------------------------------------------------- ###

### --------- Clean Results (after handling!) ---------- ###
if [[ " ${actions[@]} " =~ " delete " ]]; then
  rm -rf ~/CLEO/build_fromfile/bin ~/CLEO/build_fromfile/*.dat ~/CLEO/build_fromfile/diffs
  mkdir ~/CLEO/build_fromfile/bin
fi
### ---------------------------------------------------- ###
