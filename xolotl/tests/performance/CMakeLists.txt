#Set the package name
SET(PACKAGE_NAME "xolotl.tests.performance")

#Set the description
SET(PACKAGE_DESCRIPTION "Tests for the Xolotl performance package")

#Include directories from the source and boost binaries
include_directories(${CMAKE_SOURCE_DIR}/xolotlPerf)
include_directories(${CMAKE_SOURCE_DIR}/xolotlPerf/dummy)
include_directories(${CMAKE_SOURCE_DIR}/xolotlPerf/standard)
include_directories(${Boost_BINARY_DIRS})
include_directories(${GPTL_INCLUDE_DIRS})
include_directories(${PAPI_INCLUDE_DIRS})

#Include directories from 
if(GPTL_FOUND AND PAPI_FOUND)
    include_directories(${CMAKE_SOURCE_DIR}/xolotlPerf/standard)
    include_directories(${GPTL_INCLUDE_DIRS})
    include_directories(${PAPI_INCLUDE_DIRS})
endif(GPTL_FOUND AND PAPI_FOUND)

#Link the directories
link_directories(${CMAKE_SOURCE_DIR}/src/xolotlPerf)
link_directories(${CMAKE_SOURCE_DIR}/src/xolotlPerf/dummy)
link_directories(${CMAKE_SOURCE_DIR}/src/xolotlPerf/standard)


#if(GPTL_FOUND AND PAPI_FOUND)
#    link_directories(${CMAKE_SOURCE_DIR}/src/xolotlPerf/standard)
#endif(GPTL_FOUND AND PAPI_FOUND)

#Get the test files
file(GLOB tests *Tester.cpp)

#If boost was found, create tests
if(Boost_FOUND)
    #Make executables and link libraries for testers
    foreach(test ${tests})
        message(STATUS "Making test ${test}")
        get_filename_component(testName ${test} NAME_WE)
        add_executable(${testName} ${test})
        target_link_libraries(${testName} xolotlPerf ${GPTL_LIBRARIES} ${PAPI_LIBRARIES})
        #if(GPTL_FOUND AND PAPI_FOUND)
        #    target_link_libraries(${LIBRARY_NAME} ${GPTL_LIBRARIES} ${PAPI_LIBRARIES})
        #endif(GPTL_FOUND AND PAPI_FOUND)
        add_test(${testName} ${testName}) 
        #add a label so the tests can be run separately
        set_property(TEST ${testName} PROPERTY LABELS ${PACKAGE_NAME})   
    endforeach(test ${tests})
endif(Boost_FOUND)

