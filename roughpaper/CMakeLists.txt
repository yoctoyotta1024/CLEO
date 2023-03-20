# run cmake with " cmake -S [source dir of targerts] -B [build dir for makefile] "
# and then run cmake with this compiler:
# CXX=/opt/homebrew/bin/g++-12 cmake -S ./ -B ./build

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
set(CLEODIR "/Users/yoctoyotta1024/Documents/b1_springsummer2023/CLEO")
target_include_directories(${PROJECT_NAME} PRIVATE "${CLEODIR}/src/runmodel")
target_include_directories(${PROJECT_NAME} PRIVATE "${CLEODIR}/libs")
target_include_directories(${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR})

add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../libs" ${CMAKE_CURRENT_BINARY_DIR}/libs)
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/../src/runmodel" ${CMAKE_CURRENT_BINARY_DIR}/src/runmodel)
target_link_libraries("${PROJECT_NAME}" PRIVATE superdrop_solver initialisation runmodel)

set_target_properties(${PROJECT_NAME} PROPERTIES
  CMAKE_CXX_STANDARD_REQUIRED ON
  CMAKE_CXX_EXTENSIONS OFF
  CXX_STANDARD 20)

target_compile_options(${PROJECT_NAME} PRIVATE -Wall -g)
#target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Werror -pedantic -g)