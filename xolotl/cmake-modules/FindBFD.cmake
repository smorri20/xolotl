# Find BFD library if available.
#
# Based on FindEAVL.cmake from the XOLOTL package.
#
# Use it with
#    find_package(BFD)
#
# On some systems, libbfd won't be installed in a system location.  You can 
# specify where to search for the libbfd library using the following 
# CMake variable:
#
# BINUTILS_PREFIX    The root of the Binutils installation containing 
#                    the BFD library.
#
# Variables defined by this module:
#    
#    BFD_FOUND         System has BFD library.
#    BFD_INCLUDE_DIRS  Include directories for using the libbfd.
#    BFD_LIBRARIES     The BFD library and any needed dependencies.


include(LibFindMacros)

# Dependencies
libfind_package(BFD Libiberty)

# Search for the BFD header
find_path(BFD_INCLUDE_DIR bfd.h HINTS ${BINUTILS_PREFIX}/include)

find_library(BFD_LIBRARY bfd HINTS ${BINUTILS_PREFIX}/lib64 ${BINUTILS_PREFIX}/lib)

set(BFD_INCLUDE_DIRS ${BFD_INCLUDE_DIR})
set(BFD_LIBRARIES ${BFD_LIBRARY} ${Libiberty_LIBRARIES})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(BFD DEFAULT_MSG BFD_INCLUDE_DIR BFD_LIBRARY)

mark_as_advanced(BFD_INCLUDE_DIR BFD_LIBRARY)

