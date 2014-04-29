# This CMake file was taken, then modified, from:
#
# Try to find HDF5 headers and libraries.
#
# Usage of this module as follows:
#
#     find_package(HDF5)
#
# Variables used by this module, they can change the default behaviour and need
# to be set before calling find_package:
#
#  HDF5_PREFIX         Set this variable to the root installation of
#                      libHDF5 if the module has problems finding the
#                      proper installation path.
#
# Variables defined by this module:
#
#  HDF5_FOUND              System has HDF5 libraries and headers
#  HDF5_LIBRARIES          The HDF5 library
#  HDF5_INCLUDE_DIRS       The location of HDF5 headers


find_path(HDF5_PREFIX include/hdf5.h
    NAMES HINTS ENV HDF5_PREFIX
)

find_library(HDF5_LIBRARIES
    NAMES libhdf5.so
    HINTS ${HDF5_PREFIX}/lib
)

find_path(HDF5_INCLUDE_DIRS
    NAMES hdf5.h
    HINTS ${HDF5_PREFIX}/include
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(HDF5 DEFAULT_MSG
    HDF5_LIBRARIES
    HDF5_INCLUDE_DIRS 
)

mark_as_advanced(
    HDF5_PREFIX_DIRS
    HDF5_LIBRARIES
    HDF5_INCLUDE_DIRS 
)
