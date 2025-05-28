# module load gcc/11.2.0-gcc-11.2.0 openmpi/4.1.2-gcc-11.2.0
# CXX=mpic++ CC=mpicc cmake -S ./ -B ./build -DKokkos_ARCH_NATIVE=ON -DKokkos_ENABLE_SERIAL=ON

# set cmake version
if(NOT DEFINED CMAKE_MINIMUM_REQUIRED_VERSION)
  cmake_minimum_required(VERSION 3.18.0)
endif()

# set project name and print directory of this CMakeLists.txt (source directory of project)
project("roughpaper_scratch_test")
message(STATUS "CLEO including ${PROJECT_NAME} with PROJECT_SOURCE_DIR: ${PROJECT_SOURCE_DIR}")

#set compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -g -gdwarf-4")

# Set libraries from CLEO to link with executable
set(CLEOLIBS configuration coupldyn_fromfile cartesiandomain gridboxes initialise runcleo superdrops zarr observers)

# create primary executable "test"
add_executable(test main.cpp)

# Add directories and link libraries to target
target_link_libraries(test PRIVATE ${CLEOLIBS})
target_link_libraries(test PUBLIC Kokkos::kokkos)
target_include_directories(test PRIVATE "${CLEO_SOURCE_DIR}/libs") # CLEO libs directory
target_include_directories(test PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")

# set C++ properties for target
set_target_properties(test PROPERTIES
  CMAKE_CXX_STANDARD_REQUIRED ON
  CMAKE_CXX_EXTENSIONS ON
  CXX_STANDARD 20)
