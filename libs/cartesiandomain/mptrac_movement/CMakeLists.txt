# set cmake version
if(NOT DEFINED CMAKE_MINIMUM_REQUIRED_VERSION)
  cmake_minimum_required(VERSION 3.18.0)
endif()

# set library name and print directory of this CMakeLists.txt (source directory of library)
set(LIBNAME "mptrac_movement")
message("${LIBNAME} LIBRARY_SOURCE_DIR: ${CMAKE_CURRENT_SOURCE_DIR}")

# explicitly set library executables path to /lib in top level of build tree
set(LIB_BINARY_DIR ${CMAKE_BINARY_DIR}/lib)

### --- Create MPTRAC Library --- ###
# fetch and make MPTRAC library available
include(FetchContent)
FetchContent_Declare(
  mptrac
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  URL https://github.com/slcs-jsc/mptrac/archive/refs/heads/ww_dd.zip
  # GIT_TAG 0bff3c8e82c4b6b9b48e6472898682d1d50e7f2b
  SOURCE_SUBDIR  src
)
FetchContent_MakeAvailable(mptrac)

message(STATUS "mptrac installation in: ${CMAKE_CURRENT_BINARY_DIR}")

# file(GLOB MPTRAC_DEPENDENCY_SOURCES
#      ${mptrac_SOURCE_DIR}/Makefile
#      ${mptrac_SOURCE_DIR}/src/*.c
#      ${mptrac_SOURCE_DIR}/src/*.f90)

# add_custom_target(
#     build_mptrac_dependency
#     COMMAND make lib MPI=1 STATIC=0
#     WORKING_DIRECTORY ${mptrac_SOURCE_DIR}/src
#     DEPENDS ${MPTRAC_DEPENDENCY_SOURCES}
# )

# Create an imported target for mptrac
# # add_library(mptrac::mptrac SHARED IMPORTED)
# set_property(TARGET mptrac::mptrac PROPERTY
#     IMPORTED_LOCATION ${mptrac_SOURCE_DIR}/libmptrac.so
# )
### --- Finish creating MPTRAC Library --- ###

# Add executables and create library target
set(SOURCES
  "mptrac_transport_across_domain.cpp"
)
# must use STATIC (not(!) SHARED) lib for linking to executable if build is CUDA enabled with Kokkos
add_library("${LIBNAME}" STATIC ${SOURCES})
# add_dependencies("${LIBNAME}" build_mptrac_dependency)

# Add directories for target library
# target_include_directories(${LIBNAME} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
target_include_directories(${LIBNAME} PRIVATE "${CLEO_SOURCE_DIR}/libs") # CLEO libs directory
# target_include_directories(${LIBNAME} PRIVATE "${mptrac_SOURCE_DIR}/src") # CLEO libs directory

# Add directories and link libraries for target library
# target_include_directories(${LIBNAME} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
set(LINKLIBS superdrops gridboxes)
target_link_libraries(${LIBNAME} PUBLIC "${LINKLIBS}")
target_link_libraries(${LIBNAME} PRIVATE ${mptrac_SOURCE_DIR}/src/libmptrac.so netcdf gsl gslcblas)
target_link_libraries(${LIBNAME} PUBLIC Kokkos::kokkos)

# set specific C++ compiler options for target (optional)
#target_compile_options(${LIBNAME} PRIVATE)

# set C++ properties for target
set_target_properties(${LIBNAME} PROPERTIES
  LIBRARY_OUTPUT_DIRECTORY ${LIB_BINARY_DIR}
  CMAKE_CXX_STANDARD_REQUIRED ON
  CMAKE_CXX_EXTENSIONS ON
  CXX_STANDARD 20)
