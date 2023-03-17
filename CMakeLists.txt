#run cmake with these compiler options:
# CXX=/opt/homebrew/bin/g++-12 cmake -S ./ -B ./build

cmake_minimum_required(VERSION 3.18.0)

project(project
  LANGUAGES CXX
  DESCRIPTION "this is CLEO by Clara Bayley"
  VERSION 1.0
)

add_subdirectory(src)