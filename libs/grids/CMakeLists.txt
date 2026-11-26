# set cmake version
if(NOT DEFINED CMAKE_MINIMUM_REQUIRED_VERSION)
  cmake_minimum_required(VERSION 3.18.0)
endif()

# ${CLEO_SOURCE_DIR}/libs/grids is a shared include path for every
# grid library. Allows own nested libraries in each grid's library to
# to resolve bare includes (e.g. "cartesiandomain/domainboundaries.hpp")
target_include_directories(cleo_libs_include INTERFACE "${CLEO_SOURCE_DIR}/libs/grids")

# CLEO's Domain Options
set(valid_domain "cartesian")
if("${CLEO_DOMAIN}" IN_LIST valid_domain)
  if("${CLEO_DOMAIN}" STREQUAL "cartesian")
    add_subdirectory(cartesiandomain EXCLUDE_FROM_ALL)
  endif()
elseif("${CLEO_DOMAIN}" STREQUAL "" OR "${CLEO_DOMAIN}" STREQUAL "all") # default and all
  add_subdirectory(cartesiandomain EXCLUDE_FROM_ALL)
else()
  set(errmsg "Invalid CLEO_DOMAIN '${CLEO_DOMAIN}', please choose from: ${valid_domain}")
  message(FATAL_ERROR "${errmsg}")
endif()
