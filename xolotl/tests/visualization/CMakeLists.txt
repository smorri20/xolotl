#Set the package name
SET(PACKAGE_NAME "xolotl.tests.viz")

#Set the description
SET(PACKAGE_DESCRIPTION "Tests for the Xolotl visualization package")

#Include directories from the source and boost binaries
include_directories(${CMAKE_SOURCE_DIR}/xolotlViz)
include_directories(${CMAKE_SOURCE_DIR}/xolotlViz/dataprovider)
include_directories(${CMAKE_SOURCE_DIR}/xolotlViz/labelprovider)
include_directories(${CMAKE_SOURCE_DIR}/xolotlViz/plot)
include_directories(${Boost_BINARY_DIRS})

#Link the directories
link_directories(${CMAKE_SOURCE_DIR}/src/xolotlViz)
link_directories(${CMAKE_SOURCE_DIR}/src/xolotlViz/dataprovider)
link_directories(${CMAKE_SOURCE_DIR}/src/xolotlViz/labelprovider)
link_directories(${CMAKE_SOURCE_DIR}/src/xolotlViz/plot)

#Get the test files
file(GLOB tests *.cpp)

#If boost was found, create tests
if(Boost_FOUND)
    #Make executables and link libraries for testers
    foreach(test ${tests})
        message(STATUS "Making test ${test}")
        get_filename_component(testName ${test} NAME_WE)
        add_executable(${testName} ${test})
        target_link_libraries(${testName} xolotlData xolotlPlot)
        add_test(${testName} ${testName}) 
        #add a label so the tests can be run separately
        set_property(TEST ${testName} PROPERTY LABELS ${PACKAGE_NAME})   
    endforeach(test ${tests})
endif(Boost_FOUND)

