# run cmake with " cmake -S [source dir of targerts] -B [build dir for makefile] "
#run cmake with these compiler options:
# CXX=g++-13 cmake -S ./ -B ./build
# or CXX=g++ cmake -S ./ -B ./build

cmake_minimum_required(VERSION 3.18.0)

project(project
  LANGUAGES CXX
  DESCRIPTION "this is CLEO by Clara Bayley"
  VERSION 1.0
)

add_subdirectory(libs)
add_subdirectory(src)