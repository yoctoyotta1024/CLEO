# run cmake with " cmake -S [source dir of targerts] -B [build dir for makefile] "
# run cmake with these compiler options:
# CXX=g++ CC=gcc cmake -S ./ -B ./build

# set cmake version
cmake_minimum_required(VERSION 3.18.0)
#cmake_minimum_required(VERSION 3.21.1) # for using Kokkos with NVC++ compiler 

# setlibrary name
set(LIBNAME "zarrstorage")

# print where exactly project source dir is for this CMakeLists.txt
message("${LIBNAME} LIBRARY_SOURCE_DIR: ${CMAKE_CURRENT_SOURCE_DIR}")

# explicitly set library executables path to /lib in top level of build tree
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)

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
"massmomentsstorage.cpp"
"thermostatestorage.cpp"
"zarrstores.cpp"
)
add_library(${LIBNAME} SHARED ${SOURCES})

# add directories to include for target
target_include_directories(${LIBNAME} PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_include_directories(${LIBNAME} PUBLIC "${CMAKE_SOURCE_DIR}/libs")

# link libraries in libs to this target
set(LINKLIBS superdrop_solver sdmgridboxes)
target_link_libraries(${LIBNAME} PRIVATE ${LINKLIBS})
target_link_libraries(${LIBNAME} PUBLIC Kokkos::kokkos)

# some compilation properties
set_target_properties(${LIBNAME} PROPERTIES
  CMAKE_CXX_STANDARD_REQUIRED ON
  CMAKE_CXX_EXTENSIONS ON
  CXX_STANDARD 20)

set(compileops -Werror -Wall -pedantic -g -gdwarf-4)
#set(compileops -Werror -Wall -pedantic -g) # cannot use -gdwarf-4 flag for nvc++ compiler 

target_compile_options(${LIBNAME} PRIVATE ${compileops})