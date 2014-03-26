#Set the package name
SET(PACKAGE_NAME "xolotl.label")
#Set the description
SET(PACKAGE_DESCRIPTION "Xolotl Label Provider")
#Set the library name
SET(LIBRARY_NAME "xolotlLabel")

#Collect all header filenames in this project 
#and glob them in HEADERS
file(GLOB HEADERS *.h)

#Install the xolotl header files
install(FILES ${HEADERS} DESTINATION include)
