#!/bin/bash

# -----------------------------------------------------------------------------
# install_kokkos_tools
#
# Installs Kokkos Tools.
#
# Assumes:
#   - compiler environment has already been configured
#   - mpic++ is available in PATH
# -----------------------------------------------------------------------------

install_kokkos_tools() {

    local path2toolsrepo="$1"
    local root4tools="$2"

    if [[ -z "${path2toolsrepo}" || -z "${root4tools}" ]]; then
        echo "Usage: install_kokkos_tools <path2toolsrepo> <install-prefix>"
        return 1
    fi

    local CXX
    CXX=$(command -v mpic++)

    if [[ -z "${CXX}" ]]; then
        echo "ERROR: mpic++ not found."
        return 1
    fi

    cd "${path2toolsrepo}"

    if [[ ! -d "kokkos-tools" ]]; then
        echo "ERROR: kokkos-tools repository not found."
        return 1
    fi

    mkdir -p "${root4tools}"

    cd kokkos-tools

    mkdir -p myBuild

    cmake \
        -S . \
        -B myBuild \
        -DCMAKE_CXX_COMPILER="${CXX}" \
        -DCMAKE_INSTALL_PREFIX="${root4tools}"

    cmake --build myBuild

    cmake --install myBuild

    echo
    echo "----- INSTALLATION SUMMARY -----"

    [[ -d "${root4tools}/bin" ]] && {
        echo "bin:"
        ls "${root4tools}/bin"
    }

    [[ -d "${root4tools}/lib" ]] && {
        echo "lib:"
        ls "${root4tools}/lib"
    }

    [[ -d "${root4tools}/lib64" ]] && {
        echo "lib64:"
        ls "${root4tools}/lib64"
    }

    echo "--------------------------------"
    echo "SUCCESS: Installed to ${root4tools}"
}


### ------------ Notes on using profiler for executable ------------ ###
# example for for tools installed in /path/to/tools/kokkostools/ on macOS:
# A) see tool libraries installed in /path/to/tools/kokkostools/lib/
# B) export required tool library, e.g.
#     e.g. export KOKKOS_TOOLS_LIBS=/path/to/tools/kokkostools/lib/libkp_kernel_timer.dylib
#      or  export KOKKOS_TOOLS_LIBS=/path/to/tools/kokkostools/lib/libkp_space_time_stack.dylib
# C) run executable ./[exec].exe (kokkos initialise loads dynamic library pointers)
# D) read *.dat output
#     e.g. with kp reader
#          export DYLD_LIBRARY_PATH=/path/to/tools/kokkostools/lib/:$LD_LIBRARY_PATH
#          /path/to/tools/kokkostools/bin/kp_reader *.dat > ./bin/kp_kernel_timer.txt
#    or pipe kp_space_time_stack output durign runtime: ./[exec].exe > runtime_output.txt
# E) Also note useful debugging tool to find where program crashed (e.g. inside kernel):
#   export KOKKOS_TOOLS_LIBS=/path/to/tools/kokkostools/lib/libkp_kernel_logger.dylib
### ---------------------------------------------------------------- ###
