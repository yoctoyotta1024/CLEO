# set cmake version
if(NOT DEFINED CMAKE_MINIMUM_REQUIRED_VERSION)
  cmake_minimum_required(VERSION 3.18.0)
  #cmake_minimum_required(VERSION 3.21.1) # if using Kokkos c++ with NVC++ compiler 
endif()

# set library name
set(LIBNAME "zarr")

# print where exactly project source dir is for this CMakeLists.txt
message("${LIBNAME} LIBRARY_SOURCE_DIR: ${CMAKE_CURRENT_SOURCE_DIR}")

# explicitly set library executables path to /lib in top level of build tree
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)

# add executables to library
set(SOURCES
"fsstore.cpp"
)
add_library(${LIBNAME} SHARED ${SOURCES})

# add directories to include for target
target_include_directories(${LIBNAME} PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_include_directories(${LIBNAME} PUBLIC "${CMAKE_SOURCE_DIR}/libs")

# some compilation properties
set_target_properties(${LIBNAME} PROPERTIES
  CMAKE_CXX_STANDARD_REQUIRED ON
  CMAKE_CXX_EXTENSIONS ON
  CXX_STANDARD 20)

# set compiler options
if(NOT DEFINED compileops)
  set(compileops ${defaultcompileops})
endif()

target_compile_options(${LIBNAME} PRIVATE ${compileops})