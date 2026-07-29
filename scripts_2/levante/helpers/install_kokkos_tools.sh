#!/bin/bash

### ------------------------------------------------------- ###
### Running script successfully installs kokkos tools for gcc or intel
### compiler with openmpi on Levante, assuming you have already performed
### ```git clone git@github.com:kokkos/kokkos-tools.git```
### in ${path2toolsrepo}
### ------------------------------------------------------- ###

set -e
source /etc/profile
module purge
spack unload --all

LEVANTE_BASH_SRC="${CLEO_PATH2CLEO}/scripts_2/levante/bash/src"

path2toolsrepo=$1                       # absolute path containing the kokkos-tools repo clone
compilername=${2:-gcc}                  # compile tools with "gcc" or "intel"
root4tools=${3:-${CLEO_KOKKOSTOOLS}}    # absolute path for kokkos-tools installation

source ${LEVANTE_BASH_SRC}/levante_packages.sh

if [ "${compilername}" == "" ]
then
  echo "Bad input, please specify compiler name to build kokkos-tools with"
  exit 1
elif [ "${compilername}" == "gcc" ]
then
  module load ${levante_gcc} ${levante_gcc_openmpi}
  spack load ${levante_gcc_cmake}
elif [ "${compilername}" == "intel" ]
then
  module load ${levante_intel} ${levante_intel_openmpi}
  spack load ${levante_intel_cmake}
else
  echo "Bad input, unrecognised compiler name '${compilername}'. Must be 'gcc' or 'intel'"
  exit 1
fi

if [[ "${path2toolsrepo}" == "" || "${root4tools}" == "" ]]
then
  echo "Bad inputs, please specify absolute path for kokkos-tools repo and where you want to install the tools"
  exit 1
fi

CXX="$(command -v mpic++)"

cd ${path2toolsrepo}
if [ ! -d "./kokkos-tools" ]
then
  echo "ERROR: kokkos-tools source directory not found. Please clone it into path2toolsrepo" >&2
  exit 1
fi

mkdir -p ${root4tools}
if [ ! -d ${root4tools} ]
then
  echo "ERROR: kokkos-tools build directory not found, please make sure it exists" >&2
  exit 1
fi

cd ./kokkos-tools && mkdir ./myBuild
cmake -DCMAKE_CXX_COMPILER=${CXX} \
    -S ./ -B ./myBuild \
    -DCMAKE_INSTALL_PREFIX=${root4tools}
cd ./myBuild && make && make install

echo "----- INSTALLATION SUMMARY -----"
echo "in ${root4tools}/bin:" && ls ${root4tools}/bin
if [ -d "${root4tools}/lib" ]
then
  echo "in ${root4tools}/lib:" && ls ${root4tools}/lib
fi
if [ -d "${root4tools}/lib64" ]
then
  echo "in ${root4tools}/lib64:" && ls ${root4tools}/lib64
fi
echo "--------------------------------"
echo "SUCCESS: kokkos-tools installed in ${root4tools}"

### ------------ Notes on using profiler for executable ------------ ###
# example for tools installed in /path/to/tools/kokkostools/ on linux:
# A) see tool libraries installed in /path/to/tools/kokkostools/lib64/
# B) export required tool library, e.g.
#     e.g. export KOKKOS_TOOLS_LIBS=/path/to/tools/kokkostools/lib64/libkp_kernel_timer.so
#      or  export KOKKOS_TOOLS_LIBS=/path/to/tools/kokkostools/lib64/libkp_space_time_stack.so
# C) run executable ./[exec].exe (kokkos initialise loads dynamic library pointers)
# D) read *.dat output
#     e.g. with kp reader
#          export LD_LIBRARY_PATH=/path/to/tools/kokkostools/lib64/:$LD_LIBRARY_PATH
#          /path/to/tools/kokkostools/bin/kp_reader *.dat > ./bin/kp_kernel_timer.txt
#    or pipe kp_space_time_stack output during runtime: ./[exec].exe > runtime_output.txt
# E) Also useful debugging tool to find where program crashed (e.g. inside kernel):
#   export KOKKOS_TOOLS_LIBS=/path/to/tools/kokkostools/lib64/libkp_kernel_logger.so
### ---------------------------------------------------------------- ###
