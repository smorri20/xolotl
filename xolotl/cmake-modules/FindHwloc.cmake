# Try to find Hwloc headers and libraries.
#
# Usage of this module as follows:
#
#     find_package(Hwloc [REQUIRED])
#
# Variables used by this module, they can change the default behaviour and need
# to be set before calling find_package:
#
#  HWLOC_ROOT         Set this variable to the root of a Hwloc installation
#                      if the module has problems finding the
#                      proper installation path.
#
# Variables defined by this module:
#
#  Hwloc_FOUND            The system has Hwloc libraries and headers
#  Hwloc_LIBRARIES        Specification of Hwloc libraries to link against.
#  Hwloc_INCLUDE_DIRS     The location of Hwloc headers


find_path(HWLOC_ROOT include/hwloc.h
    NAMES HINTS ENV HWLOC_ROOT
)

find_library(Hwloc_LIBRARIES
    NAMES libhwloc.a libhwloc.so libhwloc.dylib
    HINTS ${HWLOC_ROOT}/lib
)

find_path(Hwloc_INCLUDE_DIRS
    NAMES include/hwloc.h
    HINTS ${HWLOC_ROOT}/include
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Hwloc DEFAULT_MSG
    Hwloc_LIBRARIES
    Hwloc_INCLUDE_DIRS 
)

