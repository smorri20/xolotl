# Find libunwind library if available.
#
# Based on FindEAVL.cmake from the XOLOTL package.
#
# Use it with
#    find_package(Unwind)
#
# On some systems, libunwind won't be installed in a system location.
# You can specify where to search for libunwind using the following
# CMake variable:
#
# UNWIND_PREFIX    The root of the libunwind installation.
#
# Variables defined by this module:
#    
#    Unwind_FOUND          System has libunwind library.
#    Unwind_INCLUDE_DIRS   Include directories for using the libunwind library.
#    Unwind_LIBRARIES      The libunwind library and any needed dependencies.

# Search for the libunwind header
find_path(Unwind_INCLUDE_DIR libunwind.h HINTS ${UNWIND_PREFIX}/include)

find_library(Unwind_LIBRARY unwind HINTS ${UNWIND_PREFIX}/lib)

set(Unwind_INCLUDE_DIRS ${Unwind_INCLUDE_DIR})
set(Unwind_LIBRARIES ${Unwind_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Unwind DEFAULT_MSG Unwind_INCLUDE_DIR Unwind_LIBRARY)

mark_as_advanced(Unwind_INCLUDE_DIR Unwind_LIBRARY)

