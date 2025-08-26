#!/bin/bash
#SBATCH --job-name=kokkostools
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --gpus=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=256
#SBATCH --mem=10G
#SBATCH --time=00:30:00
#SBATCH --mail-user=clara.bayley@mpimet.mpg.de
#SBATCH --mail-type=FAIL
#SBATCH --account=bm1183
#SBATCH --output=./kokkostools_out.%j.out
#SBATCH --error=./kokkostools_err.%j.out

### ------------------ Input Parameters ---------------- ###
### ------ You MUST edit these lines to set your ------- ###
### ---- build type, directories, the executable(s) ---- ###
### -------- to compile, and your python script -------- ###
### ---------------------------------------------------- ###
do_build="true"
compilername="gcc"  # must be gcc for buildtype=cuda
path2CLEO=${HOME}/CLEO/
path2build_parent=${HOME}/CLEO/build_spdtest/
build_flags="-DCLEO_COUPLED_DYNAMICS=fromfile -DCLEO_DOMAIN=cartesian \
  -DCLEO_NO_ROUGHPAPER=true -DCLEO_NO_PYBINDINGS=true"
path2kokkostools=/work/bm1183/m300950/kokkos_tools_lib/lib64/
executables="spdtest"

pythonscript=${path2CLEO}/examples/kokkostools/kokkostools.py
src_config_filename=${path2CLEO}/examples/kokkostools/src/config/kokkostools_config.yaml
postproc_filedirectory=${path2build_parent}/bin
### ---------------------------------------------------- ###
### ---------------------------------------------------- ###
### ---------------------------------------------------- ###

# ensure this directory exists
mkdir -p ${path2build_parent}

### ---- run test for different types of parallelism ---- ###
buildtypes=("cuda" "openmp" "threads" "serial")
for buildtype in "${buildtypes[@]}"
do
  path2build="${path2build_parent}/${buildtype}"
  postproc_filename="${postproc_filedirectory}/${executables}_${buildtype}.txt"
  script_args="${path2kokkostools} ${src_config_filename} ${postproc_filename} \
  --nruns=2 --do_inputfiles --do_run_executable --do_plot_results"

  echo "build type: ${buildtype}"
  echo "path2build: ${path2build}"

  ### ---------- build, compile and run example ---------- ###
  ${path2CLEO}/scripts/levante/examples/build_compile_run_plot.sh ${do_build} \
    ${buildtype} ${compilername} ${path2CLEO} ${path2build} "${build_flags}" \
    "${executables}" ${pythonscript} "${script_args}"
  ### ---------------------------------------------------- ###
done
### ---------------------------------------------------- ###
