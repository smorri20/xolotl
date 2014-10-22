#Set the package name
SET(PACKAGE_NAME "xolotl.factory")
#Set the description
SET(PACKAGE_DESCRIPTION "Xolotl Factory")
#Set the library name
SET(LIBRARY_NAME "xolotlFactory")

#Collect all header filenames in this project 
#and glob them in HEADERS
file(GLOB HEADERS *.h)

#Grab all of the source files
file(GLOB SRC *.cpp)

#Include directories
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/factory)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/advection)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/flux)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/commandline)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/reactants)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/reactants/psiclusters)
include_directories(${CMAKE_SOURCE_DIR}/xolotlPerf)

#Add a library to hold the factory code
add_library(${LIBRARY_NAME} STATIC ${SRC})
target_link_libraries(${LIBRARY_NAME})

#Install the xolotl header files
install(FILES ${HEADERS} DESTINATION include)
