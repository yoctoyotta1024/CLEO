#!/bin/bash

builddir=${1-/p/project/exaww/bayley1/CLEO/build_fromfile} # path into build directory

rm -rf ${builddir}/_deps/mptrac-*
rm ${builddir}/examples/fromfile/src/fromfile
