#Set the package name
SET(PACKAGE_NAME "xolotl.bubbles")
#Set the description
SET(PACKAGE_DESCRIPTION "Xolotl Bubbles")
#Set the library name
SET(LIBRARY_NAME "xolotlBubbles")

#Collect all header filenames in this project 
#and glob them in HEADERS
file(GLOB HEADERS *.h)

#Collect all of the cpp files in this folder 
#and the psiclusters folder and glob them in SRC
file(GLOB SRC *.cpp)

#Include this directory so that it can be built
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore)
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore/bubbles)

#Add the library
add_library(${LIBRARY_NAME} STATIC ${SRC})
target_link_libraries(${LIBRARY_NAME})

#Install the xolotl header files
install(FILES ${HEADERS} DESTINATION include)
install(TARGETS ${LIBRARY_NAME} DESTINATION lib)