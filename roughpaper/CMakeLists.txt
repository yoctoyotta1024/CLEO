# run cmake with " cmake -S [source dir of targerts] -B [build dir for makefile] "
#run cmake with these compiler options:
# CXX=/opt/homebrew/bin/g++-13 CC=/opt/homebrew/bin/gcc-13 cmake -S ./ -B ./build
# or CXX=g++ CC=gcc cmake -S ./ -B ./build

# set cmake version
cmake_minimum_required(VERSION "3.18.0")

# set project name and source directory of .cpp files
project("testing")

# print where exactly project source dir is for this CMakeLists.txt
message("PROJECT_SOURCE_DIR: ${PROJECT_SOURCE_DIR}")

# add executables to project (only one main() function!)
set(SOURCES "testing.cpp")
add_executable("${PROJECT_NAME}" ${SOURCES})

# # add directories to include for targets to us
# set(CLEODIR "/Users/yoctoyotta1024/Documents/b1_springsummer2023/CLEO")
set(CLEODIR "../")
include_directories("${CLEODIR}/libs")
# include_directories(${PROJECT_SOURCE_DIR})
# include_directories(${CLEODIR})

add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../libs/initialisation" ${CMAKE_CURRENT_BINARY_DIR}/libs/initialisation)
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../libs/superdrop_solver" ${CMAKE_CURRENT_BINARY_DIR}/libs/superdrop_solver)
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../libs/observers" ${CMAKE_CURRENT_BINARY_DIR}/libs/observers)
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../libs/sdmgridboxes" ${CMAKE_CURRENT_BINARY_DIR}/libs/sdmgridboxes)

target_link_libraries("${PROJECT_NAME}" PUBLIC superdrop_solver initialisation sdmgridboxes observers)

# Kokkos library stuff
include(FetchContent)
FetchContent_Declare(
  Kokkos
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  URL https://github.com/kokkos/kokkos/archive/refs/tags/4.0.01.tar.gz
  GIT_TAG 62d2b6c879b74b6ae7bd06eb3e5e80139c4708e6
)
FetchContent_MakeAvailable(Kokkos)

# link kokkos library to target
target_link_libraries(${PROJECT_NAME} PUBLIC Kokkos::kokkos)


set_target_properties(${PROJECT_NAME} PROPERTIES
  CMAKE_CXX_STANDARD_REQUIRED ON
  CMAKE_CXX_EXTENSIONS ON
  CXX_STANDARD 20)

target_compile_options(${PROJECT_NAME} PUBLIC -Wall -pedantic -g -gdwarf-4)
#target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Werror -pedantic -g)