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

# Find EAVL
FIND_PACKAGE(EAVL)

if(EAVL_FOUND)
    message(STATUS "EAVL includes = ${EAVL_INCLUDE_DIRS}")
    message(STATUS "EAVL config = ${EAVL_CONFIG_DIRS}")
    message(STATUS "EAVL libraries = ${EAVL_LIBRARIES}")
    include_directories(${EAVL_INCLUDE_DIRS})
    include_directories(${EAVL_INCLUDE_DIRS}/common)
    include_directories(${EAVL_INCLUDE_DIRS}/math)
    include_directories(${EAVL_INCLUDE_DIRS}/rendering)
    include_directories(${EAVL_INCLUDE_DIRS}/operations)
    include_directories(${EAVL_INCLUDE_DIRS}/fonts)
    include_directories(${EAVL_INCLUDE_DIRS}/importers)
    include_directories(${EAVL_CONFIG_DIRS})
endif(EAVL_FOUND)

# Find MESA
FIND_PACKAGE(MESA)

if(MESA_FOUND)
    message(STATUS "Mesa includes = ${MESA_INCLUDE_DIRS}")
    message(STATUS "Mesa libraries = ${MESA_LIBRARIES}")
    include_directories(${MESA_INCLUDE_DIRS})
endif(MESA_FOUND)

#Add the library
add_library(${LIBRARY_NAME} SHARED ${SRC})
target_link_libraries(${LIBRARY_NAME} ${EAVL_LIBRARIES})
target_link_libraries(${LIBRARY_NAME} ${MESA_LIBRARIES})

#Install the xolotl header files
install(FILES ${HEADERS} DESTINATION include)
install(TARGETS ${LIBRARY_NAME} DESTINATION lib)
