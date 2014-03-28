#Set the package name
SET(PACKAGE_NAME "xolotl.viz")
#Set the description
SET(PACKAGE_DESCRIPTION "Xolotl Visualization")
#Set the library name
SET(LIBRARY_NAME "xolotlViz")

#Collect all header filenames in this project 
#and glob them in HEADERS
file(GLOB HEADERS dataprovider/*.h labelprovider/*.h plot/*.h)

#Collect all of the cpp files in this folder 
#and the psiclusters folder and glob them in SRC
file(GLOB SRC dataprovider/*.cpp plot/*.cpp)

#Add the library
add_library(${LIBRARY_NAME} SHARED ${SRC})

#Install the xolotl header files
install(FILES ${HEADERS} DESTINATION include)
install(TARGETS ${LIBRARY_NAME} DESTINATION lib)
