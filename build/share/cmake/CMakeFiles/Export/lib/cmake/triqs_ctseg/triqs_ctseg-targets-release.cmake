#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "triqs_ctseg::triqs_ctseg_c" for configuration "Release"
set_property(TARGET triqs_ctseg::triqs_ctseg_c APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(triqs_ctseg::triqs_ctseg_c PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libtriqs_ctseg_c.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS triqs_ctseg::triqs_ctseg_c )
list(APPEND _IMPORT_CHECK_FILES_FOR_triqs_ctseg::triqs_ctseg_c "${_IMPORT_PREFIX}/lib/libtriqs_ctseg_c.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
