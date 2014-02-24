#Set the package name
SET(PACKAGE_NAME "xolotl.performance")
#Set the description
SET(PACKAGE_DESCRIPTION "Xolotl Perf")
#Set the library name
SET(LIBRARY_NAME "xolotlPerf")

#Include headers in order to build
include_directories(${CMAKE_SOURCE_DIR}/xolotlPerf)
include_directories(${CMAKE_SOURCE_DIR}/xolotlPerf/dummy)
include_directories(${CMAKE_SOURCE_DIR}/xolotlPerf/standard)

# Find GPTL
FIND_PACKAGE(GPTL)

# Find PAPI
FIND_PACKAGE(PAPI)

if(GPTL_FOUND AND PAPI_FOUND)
    message(STATUS "PAPI includes = ${PAPI_INCLUDE_DIRS}")
    message(STATUS "PAPI libraries = ${PAPI_LIBRARIES}")
    message(STATUS "GPTL includes = ${GPTL_INCLUDE_DIRS}")
    message(STATUS "GPTL libraries = ${GPTL_LIBRARIES}")

#Collect all header filenames in this project 
#and glob them in HEADERS
file(GLOB HEADERS *.h dummy/*.h standard/*.h)

#Collect all of the cpp files in this folder 
#and glob them in SRC
file(GLOB SRC *.cpp dummy/*.cpp standard/*.cpp)

#Add the library
add_library(${LIBRARY_NAME} SHARED ${SRC})
target_link_libraries(${LIBRARY_NAME} ${GPTL_LIBRARIES} ${PAPI_LIBRARIES})
else(GPTL_FOUND AND PAPI_FOUND)

#Collect all header filenames in this project 
#and glob them in HEADERS
file(GLOB HEADERS *.h dummy/*.h)

#Collect all of the cpp files in this folder 
#and glob them in SRC
file(GLOB SRC *.cpp dummy/*.cpp)

#Add the library
add_library(${LIBRARY_NAME} SHARED ${SRC})

endif(GPTL_FOUND AND PAPI_FOUND)

#Install the xolotl header files
install(FILES ${HEADERS} DESTINATION include)
install(TARGETS ${LIBRARY_NAME} DESTINATION lib)





