#Set the package name
SET(PACKAGE_NAME "xolotl.memusage")
#Set the description
SET(PACKAGE_DESCRIPTION "Xolotl Memory Usage")
#Set the library name
SET(LIBRARY_NAME "xolotlMemUsage")

#Include headers in order to build
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_BINARY_DIR})

# Always include the Dummy (stub), standard, and common classes.
file(GLOB HEADERS *.h common/*.h dummy/*.h standard/*.h profileproc/*.h profilenode/*.h)
file(GLOB SRC *.cpp common/*.cpp dummy/*.cpp standard/*.cpp profileproc/*.cpp profilenode/*.cpp)

# Specify the library to build
add_library(${LIBRARY_NAME} STATIC ${SRC})

#Install the xolotl header files
install(FILES ${HEADERS} DESTINATION include)
install(TARGETS ${LIBRARY_NAME} DESTINATION lib)

# Configure file for the memory usage library.
# Note: must do this after all the other checks, or else the 
# contents of the generated file won't take the check results into account.
configure_file( "${CMAKE_CURRENT_SOURCE_DIR}/memUsageConfig.h.in" "${CMAKE_CURRENT_BINARY_DIR}/memUsageConfig.h")


