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

#Add the dummy and standard subdirectories
#add_subdirectory(dummy)
#add_subdirectory(standard)

# Find GPTL
FIND_PACKAGE(GPTL)

# Find PAPI
FIND_PACKAGE(PAPI)

if(GPTL_FOUND AND PAPI_FOUND)
    message(STATUS "PAPI includes = ${PAPI_INCLUDE_DIRS}")
    message(STATUS "PAPI libraries = ${PAPI_LIBRARIES}")
    message(STATUS "GPTL includes = ${GPTL_INCLUDE_DIRS}")
    message(STATUS "GPTL libraries = ${GPTL_LIBRARIES}")
    
    include_directories(${GPTL_INCLUDE_DIRS})
    include_directories(${PAPI_INCLUDE_DIRS})

#Collect all header filenames in this project 
#and glob them in HEADERS
#file(GLOB HEADERS *.h standard/*.h dummy/*.h)
file(GLOB HEADERS EventCounter.h Timer.h standard/GPTLEventCounter.h standard/GPTLTimer.h dummy/DummyEventCounter.h dummy/DummyTimer.h)

#Collect all of the cpp files in this folder 
#and glob them in SRC
#file(GLOB SRC standard/*.cpp dummy/*.cpp)
file(GLOB SRC standard/GPTLEventCounter.cpp standard/GPTLTimer.cpp dummy/DummyEventCounter.cpp dummy/DummyTimer.cpp)

#Add the library
add_library(${LIBRARY_NAME} SHARED ${SRC})
target_link_libraries(${LIBRARY_NAME} xolotlPerf ${GPTL_LIBRARIES} ${PAPI_LIBRARIES})

#Install the xolotl header files
install(FILES ${HEADERS} DESTINATION include)
install(TARGETS ${LIBRARY_NAME} DESTINATION lib)

else(GPTL_FOUND AND PAPI_FOUND)

#Collect all header filenames in this project 
#and glob them in HEADERS
file(GLOB HEADERS *.h dummy/*.h)

#Collect all of the cpp files in this folder 
#and glob them in SRC
file(GLOB SRC dummy/*.cpp)

#Add the library
add_library(${LIBRARY_NAME} SHARED ${SRC})
target_link_libraries(${LIBRARY_NAME} xolotlPerf)

#Install the xolotl header files
install(FILES ${HEADERS} DESTINATION include)
install(TARGETS ${LIBRARY_NAME} DESTINATION lib)

endif(GPTL_FOUND AND PAPI_FOUND)







