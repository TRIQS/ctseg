# This file allows other CMake Projects to find us
# We provide general project information
# and reestablish the exported CMake Targets

# Multiple inclusion guard
if(NOT triqs_ctseg_FOUND)
set(triqs_ctseg_FOUND TRUE)
set_property(GLOBAL PROPERTY triqs_ctseg_FOUND TRUE)

# version
set(triqs_ctseg_VERSION 3.1.0 CACHE STRING "triqs_ctseg version")
set(triqs_ctseg_GIT_HASH 6983f08f8c4f0504fceacf77adc2c4f81fce33a4 CACHE STRING "triqs_ctseg git hash")

# Root of the installation
set(triqs_ctseg_ROOT /mnt/home/nkavokine/triqs_unst CACHE STRING "triqs_ctseg root directory")

## Find the target dependencies
#function(find_dep)
#  get_property(${ARGV0}_FOUND GLOBAL PROPERTY ${ARGV0}_FOUND)
#  if(NOT ${ARGV0}_FOUND)
#    find_package(${ARGN} REQUIRED HINTS /mnt/home/nkavokine/triqs_unst)
#  endif()
#endfunction()
#find_dep(depname 1.0)

# Include the exported targets of this project
include(/mnt/home/nkavokine/triqs_unst/lib/cmake/triqs_ctseg/triqs_ctseg-targets.cmake)

message(STATUS "Found triqs_ctseg-config.cmake with version 3.1.0, hash = 6983f08f8c4f0504fceacf77adc2c4f81fce33a4, root = /mnt/home/nkavokine/triqs_unst")

# Was the Project built with Documentation?
set(triqs_ctseg_WITH_DOCUMENTATION OFF CACHE BOOL "Was triqs_ctseg build with documentation?")

# Was the Project built with PythonSupport?
set(triqs_ctseg_WITH_PYTHON_SUPPORT OFF CACHE BOOL "Was triqs_ctseg build with python support?")
if(OFF)
  set(triqs_ctseg_MODULE_DIR /mnt/home/nkavokine/triqs_unst/lib/python3.9/site-packages CACHE BOOL "The triqs_ctseg python module directory")
endif()

endif()
