#Set the package name
SET(PACKAGE_NAME "xolotl.memusage")
#Set the description
SET(PACKAGE_DESCRIPTION "Xolotl Memory Usage")
#Set the library name
SET(LIBRARY_NAME "xolotlMemUsage")

#Include headers in order to build
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_BINARY_DIR})

# Always use the sources and headers common to all registry types.
file(GLOB COMMONHEADERS *.h)
file(GLOB COMMONSRC *.cpp)

# Always include the Dummy (stub) classes.
file(GLOB DUMMYHEADERS *.h dummy/*.h)
file(GLOB DUMMYSRC *.cpp dummy/*.cpp)

# Always include the Standard classes, since we will always build them.
set(STD_HEADERS standard/StdHandlerRegistry.h standard/MemUsageSampler.h)
set(STD_SRC standard/StdHandlerRegistry.cpp standard/StatmSampler.cpp)

set(HEADERS ${COMMONHEADERS} ${DUMMYHEADERS} ${STD_HEADERS})
set(SRC ${COMMONSRC} ${DUMMYSRC} ${STD_SRC})

message(STATUS "xolotlMemUsage: are we even getting here?")


# Specify the library to build
add_library(${LIBRARY_NAME} STATIC ${SRC})

#Install the xolotl header files
install(FILES ${HEADERS} DESTINATION include)
install(TARGETS ${LIBRARY_NAME} DESTINATION lib)


# Configure file for the memory usage library.
# Note: must do this after all the other checks, or else the 
# contents of the generated file won't take the check results into account.
configure_file( "${CMAKE_CURRENT_SOURCE_DIR}/memUsageConfig.h.in" "${CMAKE_CURRENT_BINARY_DIR}/memUsageConfig.h")


