# run cmake with " cmake -S [source dir of targerts] -B [build dir for makefile] "
# and then run cmake with this compiler:
# CXX=/opt/homebrew/bin/g++-12 cmake -S ./ -B ./build

# set cmake version
cmake_minimum_required(VERSION "3.18.0")

# set project name and source directory of .cpp files
project("thermo_solver")

# print where exactly project source dir is for this CMakeLists.txt
message("PROJECT_SOURCE_DIR: ${PROJECT_SOURCE_DIR}")

set(LIBRARY_OUTPUT_PATH  ${CMAKE_BINARY_DIR}/lib)

# add executables to project (only one main() function!)
set(SOURCES 
  "cvodethermosolver.cpp"
  "differential_functions.cpp"
  )
add_library("${PROJECT_NAME}" SHARED ${SOURCES})

# # add directories to include for targets to us
target_include_directories(${PROJECT_NAME} PRIVATE ${PROJECT_SOURCE_DIR})

set(INITLIB /Users/yoctoyotta1024/Documents/b1_springsummer2023/CLEO/libs)
target_include_directories(${PROJECT_NAME} PRIVATE ${INITLIB})

# CVODE libaries stuff
# specify libraries to link against
set(LIBDIR /usr/local/cvode-6/lib)
find_library(CVODE_LIB NAMES sundials_cvode HINTS "${LIBDIR}")
find_library(NVECS_LIB NAMES sundials_nvecserial HINTS "${LIBDIR}")

# find_library(MANYVECS_LIB NAMES sundials_nvecmanyvector  HINTS "${LIBDIR}")
set(SUNDIALS_LIBS ${CVODE_LIB} ${NVECS_LIB})
target_link_libraries(${PROJECT_NAME} PRIVATE ${SUNDIALS_LIBS})

# CVODE include libraries stuff
target_include_directories(${PROJECT_NAME} PUBLIC /usr/local/cvode-6/include)

set_target_properties(${PROJECT_NAME} PROPERTIES
  CMAKE_CXX_STANDARD_REQUIRED ON
  CMAKE_CXX_EXTENSIONS OFF
  CXX_STANDARD 20)

target_compile_options(${PROJECT_NAME} PRIVATE -Wall -g)
#target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Werror -pedantic -g)