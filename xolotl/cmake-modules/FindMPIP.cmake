# Find mpiP lightweight MPI profiling library if available.
#
# Based on FindEAVL.cmake, from the XOLOTL package.
#
# Use it with
#     find_package(MPIP)
#
# On most systems, mpiP won't be in a system library.  You can specify
# where to search for mpiP using the following CMake variables:
#
# MPIP_PREFIX       Set this variable to the root of the mpiP installation.
#
# Variables defined by this module:
#
#  MPIP_FOUND              System has mpiP libraries
#  MPIP_LIBRARIES          Linker path flags and library flags to link mpiP

include(LibFindMacros)

# Dependencies
libfind_package(MPIP BFD)
libfind_package(MPIP Unwind)

find_library(MPIP_LIBRARY
    NAMES libmpiP.a
    HINTS ${MPIP_PREFIX}/lib
)

set(MPIP_LIBRARIES ${MPIP_LIBRARY} ${BFD_LIBRARIES} ${Unwind_LIBRARIES})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MPIP DEFAULT_MSG MPIP_LIBRARY)

mark_as_advanced(MPIP_LIBRARY)

