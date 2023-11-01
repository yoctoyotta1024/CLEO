# module load gcc && CXX=g++ CC=gcc cmake -S ./ -B ./build -DKokkos_ARCH_NATIVE=ON -DKokkos_ENABLE_SERIAL=ON

# set cmake version
cmake_minimum_required(VERSION 3.18.0)
# cmake_minimum_required(VERSION 3.21.1) # if using Kokkos c++ with NVC++ compiler 

# set project name and source directory of .cpp files
project("test")

# print where exactly project source dir is for this CMakeLists.txt
message("PROJECT_SOURCE_DIR: ${PROJECT_SOURCE_DIR}")

# Kokkos library stuff
set(CMAKE_CXX_FLAGS_B4 "${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_STANDARD_B4 "${CMAKE_CXX_STANDARD}")
include(FetchContent)
FetchContent_Declare(
  Kokkos
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  URL https://github.com/kokkos/kokkos/archive/refs/tags/4.0.01.tar.gz
  GIT_TAG 62d2b6c879b74b6ae7bd06eb3e5e80139c4708e6
)
set(CMAKE_CXX_STANDARD "20")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
FetchContent_MakeAvailable(Kokkos)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_B4}")
set(CMAKE_CXX_STANDARD "${CMAKE_CXX_STANDARD_B4}")

# # Set libraries and compiler warnings to use
# set(SDMLIBS superdrops)

# default compiler options to use in subdirectories
set(defaultcompileops -Werror -Wall -pedantic -g -gdwarf-4)   # good for g++ compiler
# set(defaultcompileops -Werror -Wall -pedantic -g)           # -gdwarf-4 flag not possible for nvc++ compiler 

# primary executable for CLEO sdm
set(test "test")
add_executable(${test} "main.cpp")
# target_link_libraries(${test} PRIVATE "${SDMLIBS}")
target_link_libraries(${test} PUBLIC Kokkos::kokkos)
target_compile_options(${test} PRIVATE ${compileops})

# set compiler properties for target(s)
set_target_properties(${test} PROPERTIES
  CMAKE_CXX_STANDARD_REQUIRED ON
  CMAKE_CXX_EXTENSIONS ON
  CXX_STANDARD 20)