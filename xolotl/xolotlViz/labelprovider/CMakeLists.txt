#Set the package name
SET(PACKAGE_NAME "xolotl.label")
#Set the description
SET(PACKAGE_DESCRIPTION "Xolotl Label Provider")
#Set the library name
SET(LIBRARY_NAME "xolotlLabel")

#Collect all header filenames in this project 
#and glob them in HEADERS
file(GLOB HEADERS *.h)

#Collect all of the cpp files in this folder 
#and glob them in SRC
file(GLOB SRC *.cpp)

#Add the library
add_library(${LIBRARY_NAME} SHARED ${SRC})

#Install the xolotl header files
install(FILES ${HEADERS} DESTINATION include)
install(TARGETS ${LIBRARY_NAME} DESTINATION lib)
