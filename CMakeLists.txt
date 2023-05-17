# run cmake with " cmake -S [source dir of targerts] -B [build dir for makefile] "
#run cmake with these compiler options:
# CXX=g++-13 CC=gcc-13 cmake -S ./ -B ./build
# CXX=g++ CC=gcc cmake -S ./ -B ./build
# CXX=g++-13 CC=gcc-13 cmake -S ./ -B ./build -DKokkos_ENABLE_OPENMP=ON -DKokkos_ARCH_NATIVE=ON
# -DKokkos_ENABLE_CUDA=ON -DKokkos_ENABLE_CUDA_LAMBDA=ON

cmake_minimum_required(VERSION 3.18.0)

project(project
  LANGUAGES CXX
  DESCRIPTION "this is CLEO by Clara Bayley"
  VERSION 1.0
)

add_subdirectory(libs)
add_subdirectory(src)