#Set the package name
SET(PACKAGE_NAME "xolotl.memusage")
#Set the description
SET(PACKAGE_DESCRIPTION "Xolotl Memory Usage")
#Set the library name
SET(LIBRARY_NAME "xolotlMemUsage")

#Include headers in order to build
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_BINARY_DIR})

# Specify files we always need to build/depend on, regardless of anything
# else that is available/supported.
file(GLOB COMMON_HEADERS *.h common/*.h dummy/*.h)
file(GLOB COMMON_SRCS *.cpp common/*.cpp dummy/*.cpp)

message(STATUS "memusage COMMON_HEADERS=${COMMON_HEADERS}")
message(STATUS "memusage COMMON_SRCS=${COMMON_SRCS}")

include(CheckIncludeFiles)

# Check if we have sys/sysinfo.h.  We use this on Linux for per-node data.
check_include_files(sys/sysinfo.h HAVE_SYSINFO)

if(HAVE_SYSINFO)
    set(HAVE_PER_NODE_DATA_SOURCE 1)
    file(GLOB SYSINFO_HEADERS common/SysInfo/*.h summarynode/*.h profilenode/*.h)
    file(GLOB SYSINFO_SRCS common/SysInfo/*.cpp summarynode/*.cpp profilenode/*.cpp)
endif()

# Check if we have /proc.  We assume that if we do, we will have
# /proc/<pid>/statm for our per-process data source.
message(STATUS "Looking for /proc")
if(EXISTS "/proc")
    message(STATUS "Looking for /proc - found")
    set(HAVE_STATM TRUE)
else()
    message(STATUS "Looking for /proc - not found")
    set(HAVE_STATM FALSE)
endif()

if(HAVE_STATM)
    set(HAVE_PER_PROC_DATA_SOURCE 1)
    file(GLOB STATM_HEADERS common/Statm/*.h summaryproc/*.h profileproc/*.h)
    file(GLOB STATM_SRCS common/Statm/*.cpp summaryproc/*.cpp profileproc/*.cpp)
endif()

# Specify headers and source files to build/depend on.
set(HEADERS ${COMMON_HEADERS} ${STATM_HEADERS} ${SYSINFO_HEADERS})
set(SRC ${COMMON_SRCS} ${STATM_SRCS} ${SYSINFO_SRCS})

# Specify the library to build
add_library(${LIBRARY_NAME} STATIC ${SRC})

#Install the xolotl header files
install(FILES ${HEADERS} DESTINATION include)
install(TARGETS ${LIBRARY_NAME} DESTINATION lib)

# Configure file for the memory usage library.
# Note: must do this after all the other checks, or else the 
# contents of the generated file won't take the check results into account.
configure_file( "${CMAKE_CURRENT_SOURCE_DIR}/memUsageConfig.h.in" "${CMAKE_CURRENT_BINARY_DIR}/memUsageConfig.h")


