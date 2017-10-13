# Try to find Kokkos headers and libraries.
#
# Usage of this module as follows:
#
#     find_package(Kokkos [REQUIRED])
#
# Variables used by this module, they can change the default behaviour and need
# to be set before calling find_package:
#
#  KOKKOS_ROOT         Set this variable to the root of a Kokkos installation
#                      if the module has problems finding the
#                      proper installation path.
#
# Variables defined by this module:
#
#  Kokkos_FOUND            The system has Kokkos libraries and headers
#  Kokkos_LIBRARIES        Specification of Kokkos libraries to link against.
#  Kokkos_INCLUDE_DIRS     The location of Kokkos headers


find_path(KOKKOS_ROOT include/Kokkos_Core.hpp
    NAMES HINTS ENV KOKKOS_ROOT
)

find_library(Kokkos_LIBRARIES
    NAMES libkokkos.a
    HINTS ${KOKKOS_ROOT}/lib
)

find_path(Kokkos_INCLUDE_DIRS
    NAMES Kokkos_Core.hpp
    HINTS ${KOKKOS_ROOT}/include
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Kokkos DEFAULT_MSG
    Kokkos_LIBRARIES
    Kokkos_INCLUDE_DIRS 
)

