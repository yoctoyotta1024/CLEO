# run cmake with " cmake -S [source dir of targerts] -B [build dir for makefile] "
#run cmake with these compiler options:
# CXX=/opt/homebrew/bin/g++-12 cmake -S ./ -B ./build
# or CXX=gcc cmake -S ./ -B ./build

# set cmake version
cmake_minimum_required(VERSION "3.18.0")

# setlibrary name
set(LIBNAME "thermo_solver")

# print where exactly project source dir is for this CMakeLists.txt
message("${LIBNAME} LIBRARY_SOURCE_DIR: ${CMAKE_CURRENT_SOURCE_DIR}")

# explicitly set library executables path to /lib in top level of build tree
set(LIBRARY_OUTPUT_PATH  ${CMAKE_BINARY_DIR}/lib)

# add executables to library
set(SOURCES 
  "cvodethermosolver.cpp"
  "differential_functions.cpp"
  )
add_library("${LIBNAME}" SHARED ${SOURCES})

# add directories to include for targets
target_include_directories(${LIBNAME} PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_include_directories(${LIBNAME} PUBLIC "${CMAKE_SOURCE_DIR}/libs")

# CVODE libaries stuff
include(FetchContent)
FetchContent_Declare(
  cvodes
  URL https://github.com/LLNL/sundials/releases/download/v6.5.0/cvodes-6.5.0.tar.gz
  URL_HASH MD5=dc0cd29340599f0f8b027c6f1e0107ba
)
FetchContent_MakeAvailable(cvodes)

# link cvode library to target
target_link_libraries(${LIBNAME} PUBLIC SUNDIALS::cvodes_static)

# some compilation properties
set_target_properties(${LIBNAME} PROPERTIES
  CMAKE_CXX_STANDARD_REQUIRED ON
  CMAKE_CXX_EXTENSIONS ON
  CXX_STANDARD 20)

target_compile_options(${LIBNAME} PRIVATE -Wall -g)
#target_compile_options(${LIBNAME} PRIVATE -Wall -Werror -pedantic -g)