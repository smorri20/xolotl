#Set the package name
SET(PACKAGE_NAME "xolotl.performance")
#Set the description
SET(PACKAGE_DESCRIPTION "Xolotl Perf")
#Set the library name
SET(LIBRARY_NAME "xolotlPerf")

#Include headers in order to build
include_directories(${CMAKE_SOURCE_DIR}/xolotlCore)
include_directories(${CMAKE_SOURCE_DIR}/xolotlPerf)
include_directories(${CMAKE_SOURCE_DIR}/xolotlPerf/dummy)

# Find GPTL
FIND_PACKAGE(GPTL)

# Find PAPI
FIND_PACKAGE(PAPI)

# Always include the Dummy (stub) classes.
file(GLOB HEADERS *.h dummy/*.h)
file(GLOB SRC *.cpp dummy/*.cpp)


if(GPTL_FOUND AND PAPI_FOUND)
    message(STATUS "PAPI includes = ${PAPI_INCLUDE_DIRS}")
    message(STATUS "PAPI libraries = ${PAPI_LIBRARIES}")
    message(STATUS "GPTL includes = ${GPTL_INCLUDE_DIRS}")
    message(STATUS "GPTL libraries = ${GPTL_LIBRARIES}")
    
    include_directories(${CMAKE_SOURCE_DIR}/xolotlPerf/standard)
    include_directories(${GPTL_INCLUDE_DIRS})
    include_directories(${PAPI_INCLUDE_DIRS})

    # Add the standard (GPTL+PAPI) classes to the HEADERS and SRCS.
    file(GLOB STDHEADERS standard/*.h)
    set(HEADERS ${HEADERS} ${STDHEADERS})
    file(GLOB STDSRC standard/*.cpp)
    set(SRC ${SRC} ${STDSRC})

    # Make sure we link with GPTL and with PAPI.
    # TODO GPTL might not have been built with PAPI, so 
    # do we want to require it?
    set(MAYBE_GPTL_PAPI ${GPTL_LIBRARIES} ${PAPI_LIBRARIES})

endif(GPTL_FOUND AND PAPI_FOUND)


#Add the library
add_library(${LIBRARY_NAME} SHARED ${SRC})
target_link_libraries(${LIBRARY_NAME} ${MAYBE_GPTL_PAPI})

#Install the xolotl header files
install(FILES ${HEADERS} DESTINATION include)
install(TARGETS ${LIBRARY_NAME} DESTINATION lib)



