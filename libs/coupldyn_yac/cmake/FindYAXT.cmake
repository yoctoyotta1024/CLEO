find_path(YAXT_C_INCLUDE_DIR
  NAMES yaxt.h
  DOC "YAXT C include dir")

find_library(YAXT_C_LIBRARY
  NAMES libyaxt_c.a
  DOC "YAXT C Library")

mark_as_advanced(YAXT_C_INCLUDE_DIR
  YAXT_C_LIBRARY)

# CLEO is C / C++ only and therefore links YAXT::YAXT_C. YAXT's Fortran interface is
# optional: it is absent from C-only YAXT builds (e.g. one configured with FC=no).
find_path(YAXT_Fortran_INCLUDE_DIR
  NAMES yaxt.mod
  HINTS ${YAXT_C_INCLUDE_DIR})

find_library(YAXT_Fortran_LIBRARY
  NAMES libyaxt.a
  DOC "YAXT Fortran Library (optional)")

mark_as_advanced(YAXT_Fortran_INCLUDE_DIR
  YAXT_Fortran_LIBRARY)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(YAXT
REQUIRED_VARS YAXT_C_LIBRARY YAXT_C_INCLUDE_DIR
)

if(YAXT_FOUND)
  if(NOT TARGET YAXT::YAXT_C)
    add_library(YAXT::YAXT_C INTERFACE IMPORTED)
    set_property(TARGET YAXT::YAXT_C
      PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${YAXT_C_INCLUDE_DIR}")
    set_property(TARGET YAXT::YAXT_C
      PROPERTY INTERFACE_LINK_LIBRARIES "${YAXT_C_LIBRARY}")
  endif()

  if(YAXT_Fortran_LIBRARY AND NOT TARGET YAXT::YAXT_Fortran)
    add_library(YAXT::YAXT_Fortran INTERFACE IMPORTED)
    set_property(TARGET YAXT::YAXT_Fortran
      PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${YAXT_Fortran_INCLUDE_DIR}")
    set_property(TARGET YAXT::YAXT_Fortran
      PROPERTY INTERFACE_LINK_LIBRARIES "${YAXT_Fortran_LIBRARY}")
    target_link_libraries(YAXT::YAXT_Fortran INTERFACE YAXT::YAXT_C)
  endif()
endif()
