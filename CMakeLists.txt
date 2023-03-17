#run cmake with these compiler options:
# CXX=/opt/homebrew/bin/g++-12 cmake -S ./ -B ./build

cmake_minimum_required(VERSION 3.23.2)
project(example
  LANGUAGES CXX
  VERSION 1.0
)

add_subdirectory(src)