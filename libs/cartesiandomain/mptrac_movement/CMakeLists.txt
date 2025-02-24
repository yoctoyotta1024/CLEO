# set cmake version
if(NOT DEFINED CMAKE_MINIMUM_REQUIRED_VERSION)
  cmake_minimum_required(VERSION 3.18.0)
endif()

# set library name and print directory of this CMakeLists.txt (source directory of library)
set(LIBNAME "cartesian_mptrac_movement")
message("${LIBNAME} LIBRARY_SOURCE_DIR: ${CMAKE_CURRENT_SOURCE_DIR}")

# explicitly set library executables path to /lib in top level of build tree
set(LIB_BINARY_DIR ${CMAKE_BINARY_DIR}/lib)

### --- Create MPTRAC src C library --- ###
# fetch and make MPTRAC available
include(FetchContent)
FetchContent_Declare(
  mptrac
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  URL https://github.com/slcs-jsc/mptrac/archive/refs/heads/ww_dd.zip
  # GIT_TAG 0bff3c8e82c4b6b9b48e6472898682d1d50e7f2b
)
FetchContent_MakeAvailable(mptrac)

message(STATUS "MPTRAC installation in: ${mptrac_SOURCE_DIR}/src")

enable_language(C)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -g -DHAVE_INLINE -fno-common -fshort-enums -fopenmp")
set(mptrac_dependency_libs netcdf gsl gslcblas)

add_library(mptrac_lib SHARED "${mptrac_SOURCE_DIR}/src/mptrac.c")
target_include_directories(mptrac_lib PRIVATE "${mptrac_SOURCE_DIR}/src") # MPTRAC src
target_link_libraries(mptrac_lib PRIVATE ${mptrac_dependency_libs}) # MPTRAC dependencies

# note: -DMPI enables MPI for MPTRAC (may not be compatible with static library?)
target_compile_definitions(mptrac_lib PRIVATE
  MPI=${MPTRAC_MPI}
  VERSION="${MPTRAC_VERSION}"
)

set_target_properties(mptrac_lib PROPERTIES
  C_STANDARD 99
  C_STANDARD_REQUIRED ON
  C_EXTENSIONS OFF
)
### --------------------------------- ###

# Add executables and create library target
set(SOURCES
  "mptrac_transport_across_domain.cpp"
)
# must use STATIC (not(!) SHARED) lib for linking to executable if build is CUDA enabled with Kokkos
add_library("${LIBNAME}" STATIC ${SOURCES})
add_dependencies("${LIBNAME}" mptrac_lib) # library depend on custom target for building MPTRAC src

# Add directories for target library
# target_include_directories(${LIBNAME} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
target_include_directories(${LIBNAME} PRIVATE "${CLEO_SOURCE_DIR}/libs") # CLEO libs
target_include_directories(${LIBNAME} PUBLIC "${mptrac_SOURCE_DIR}/src") # MPTRAC src

# Link libraries for target library
set(LINKLIBS superdrops gridboxes)
target_link_libraries(${LIBNAME} PUBLIC "${LINKLIBS}")
target_link_libraries(${LIBNAME} PRIVATE mptrac_lib)
target_link_libraries(${LIBNAME} PUBLIC Kokkos::kokkos)

# set specific C++ compiler options for target (optional)
#target_compile_options(${LIBNAME} PRIVATE)

# set C++ properties for target
set_target_properties(${LIBNAME} PROPERTIES
  LIBRARY_OUTPUT_DIRECTORY ${LIB_BINARY_DIR}
  CMAKE_CXX_STANDARD_REQUIRED ON
  CMAKE_CXX_EXTENSIONS ON
  CXX_STANDARD 20)
