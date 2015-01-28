# Find Libiberty library if available.
#
# Based on FindEAVL.cmake from the XOLOTL package.
#
# Use it with
#    find_package(Libiberty)
#
# On some systems, the library won't be installed in a system location.  
# You can specify where to search for the library using the following 
# CMake variable:
#
# BINUTILS_PREFIX    The root of the Binutils installation containing 
#                    the Libiberty library.
#
# Variables defined by this module:
#    
#    Libiberty_FOUND        System has the Libiberty library.
#    Libiberty_LIBRARIES    The Libiberty library and any needed dependencies.


find_library(Libiberty_LIBRARY iberty HINTS ${BINUTILS_PREFIX}/lib64)

set(Libiberty_LIBRARIES ${Libiberty_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Libiberty DEFAULT_MSG Libiberty_LIBRARY)

mark_as_advanced(Libiberty_LIBRARY)

