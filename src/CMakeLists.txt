# set cmake version
if(NOT DEFINED CMAKE_MINIMUM_REQUIRED_VERSION)
  cmake_minimum_required(VERSION 3.18.0)
endif()

# set project name and print directory of this CMakeLists.txt (source directory of project)
project("runCLEO")
message("${PROJECT_NAME} PROJECT_SOURCE_DIR: ${PROJECT_SOURCE_DIR}")

# # Kokkos library stuff
# set(CMAKE_CXX_FLAGS_B4 "${CMAKE_CXX_FLAGS}")
# set(CMAKE_CXX_STANDARD_B4 "${CMAKE_CXX_STANDARD}")
# include(FetchContent)
# FetchContent_Declare(
#   Kokkos
#   DOWNLOAD_EXTRACT_TIMESTAMP TRUE
#   URL https://github.com/kokkos/kokkos/archive/refs/tags/4.2.00.tar.gz
#   GIT_TAG 71a9bcae52543bd065522bf3e41b5bfa467d8015
# )
# set(CMAKE_CXX_STANDARD "20")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
# FetchContent_MakeAvailable(Kokkos)
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_B4}")
# set(CMAKE_CXX_STANDARD "${CMAKE_CXX_STANDARD_B4}")

# Set libraries from CLEO to link with executable
set(CLEOLIBS gridboxes initialise observers runcleo superdrops zarr)

# create primary executable "runCLEO"
set(runCLEO "runCLEO")
add_executable(${runCLEO} "main.cpp")

# Add directories and libraries for target
target_include_directories(${runCLEO} PRIVATE "${CMAKE_SOURCE_DIR}/libs")
# target_include_directories(${runCLEO} PRIVATE ${PROJECT_SOURCE_DIR})
target_link_libraries(${runCLEO} PRIVATE coupldyn_fromfile cartesiandomain "${SDMLIBS}")
target_link_libraries(${runCLEO} PUBLIC Kokkos::kokkos)

# set C++ compiler options for target
target_compile_options(${runCLEO} PRIVATE ${CMAKE_CXX_FLAGS})

# set C++ properties for target
set_target_properties(${runCLEO} PROPERTIES
  CMAKE_CXX_STANDARD_REQUIRED ON
  CMAKE_CXX_EXTENSIONS ON
  CXX_STANDARD 20)
