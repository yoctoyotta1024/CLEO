# run cmake with " cmake -S [source dir of targerts] -B [build dir for makefile] "
#run cmake with these compiler options:
# CXX=g++-13 CC=gcc-13 cmake -S ./ -B ./build
# or CXX=g++ CC=gcc cmake -S ./ -B ./build

# set cmake version
#cmake_minimum_required(VERSION 3.18.0)
cmake_minimum_required(VERSION 3.21.1) # for using Kokkos with NVC++ compiler 

# setlibrary name
set(LIBNAME "initialisation")

# print where exactly project source dir is for this CMakeLists.txt
message("${LIBNAME} LIBRARY_SOURCE_DIR: ${CMAKE_CURRENT_SOURCE_DIR}")

# explicitly set library executables path to /lib in top level of build tree
set(LIBRARY_OUTPUT_PATH  ${CMAKE_BINARY_DIR}/lib)

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

# add executables to library
set(SOURCES 
  "config.cpp"
  "read_gbxboundaries.cpp"
  "read_initsuperdrops.cpp"
  "readbinary.cpp"
  )
add_library("${LIBNAME}" SHARED ${SOURCES})

# add directories to include for targets
target_include_directories(${LIBNAME} PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(${LIBNAME} PUBLIC Kokkos::kokkos)

# some compilation properties
set_target_properties(${LIBNAME} PROPERTIES
  CMAKE_CXX_STANDARD_REQUIRED ON
  CMAKE_CXX_EXTENSIONS ON
  CXX_STANDARD 20)

target_compile_options(${LIBNAME} PRIVATE -Wall -pedantic -g -gdwarf-4)
#target_compile_options(${LIBNAME} PRIVATE -Wall -Werror -pedantic -g -gdwarf-4)