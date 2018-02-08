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

# We need to know some of the configuration needed to support
# whichever Kokkos backend(s) we might use.
if(EXISTS ${KOKKOS_ROOT}/kokkos.cmake)
    # Older Kokkos revision (< 2.5)
    include("${KOKKOS_ROOT}/kokkos.cmake")
elseif(EXISTS ${KOKKOS_ROOT}/kokkos_generated_settings.cmake)
    # Newer Kokkos revision (2.5+)
    include("${KOKKOS_ROOT}/kokkos_generated_settings.cmake")
endif()

# * Serial, nothing else should be needed for our configuration.
if("${KOKKOS_DEVICES}" MATCHES "Serial")
    message(STATUS "Detected Kokkos support for Serial")
endif()

# For pthreads, we may need to add a pthreads library but shouldn't
# need to add anything to CXXFLAGS.
if("${KOKKOS_DEVICES}" MATCHES "Pthread")
    message(STATUS "Detected Kokkos support for Pthread")
endif()

# For OpenMP, we use CMake's support for figuring out what it
# needs to build OpenMP code rather than taking it from whatever
# configuration we used to configure Kokkos.
if("${KOKKOS_DEVICES}" MATCHES "OpenMP")
    message(STATUS "Detected Kokkos support for OpenMP.  Finding OpenMP.")
    find_package(OpenMP REQUIRED)
endif()

# For CUDA, we probably need to use whatever was used to configure Kokkos.
if("${KOKKOS_DEVICES}" MATCHES "Cuda")
    message(STATUS "Detected Kokkos support for CUDA")
    # The flags we used when building Kokkos probably have a C++
    # standard selection flag (e.g., -std=c++11).  We probably set this
    # already ourself, so to avoid compiler warnings, we strip it.
    string(REGEX REPLACE "--*std=[^ ]*" "" MY_KOKKOS_CXXFLAGS ${KOKKOS_CXXFLAGS})
    set(KOKKOS_CXXFLAGS "${MY_KOKKOS_CXXFLAGS}")

    # We assume that the user is compiling with nvcc_wrapper,
    # which adds the appropriate linker flags and libraries for us.
endif()

# For std::threads (once it is supported), nothing should be
# needed beyond whatever C++ standards flag we give.


