#Set the package name
SET(PACKAGE_NAME "xolotl.performance")
#Set the description
SET(PACKAGE_DESCRIPTION "Xolotl Perf")
#Set the library name
SET(LIBRARY_NAME "xolotlPerf")

#Collect all header filenames in this project 
#and glob them in HEADERS
file(GLOB HEADERS *.h)

# Add the dummy and standard subdirectory
add_subdirectory(dummy)
#add_subdirectory(standard)

#Include headers in order to build
include_directories(${CMAKE_SOURCE_DIR}/xolotlPerf)
include_directories(${CMAKE_SOURCE_DIR}/xolotlPerf/dummy)
#include_directories(${CMAKE_SOURCE_DIR}/xolotlPerf/standard)

#Add the library
add_library(${LIBRARY_NAME} SHARED ${SRC})

#Install the xolotl header files
install(FILES ${HEADERS} DESTINATION include)
install(TARGETS ${LIBRARY_NAME} DESTINATION lib)
