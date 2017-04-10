#Set the package name
SET(PACKAGE_NAME "xolotl.tests.memusage")

#Set the description
SET(PACKAGE_DESCRIPTION "Tests for the Xolotl memusage package")

# Our tests depend on Boost.
if(Boost_FOUND)

    # We will need the Boost headers.
    include_directories(${Boost_BINARY_DIRS})

    # Look for headers in the source tree and binary tree.
    include_directories(${CMAKE_SOURCE_DIR})
    include_directories(${CMAKE_BINARY_DIR})

    # Always build tests for the dummy classes.
    file(GLOB DUMMY_TEST_SRCS Dummy*Tester.cpp)

    # Always build the testers for the Standard classes that are always built
    set(COMMON_TEST_SRCS MemUsageSamplerTester.cpp MemUsageStdHandlerRegistryTester.cpp)

    # Make a list of all memusage infrastructure tests we will build
    set(tests ${DUMMY_TEST_SRCS} ${COMMON_TEST_SRCS})

    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/memUsageTestConfig.h.in" "${CMAKE_CURRENT_BINARY_DIR}/memUsageTestConfig.h")

    # Define executables for each of the testers.
    foreach(test ${tests})
        message(STATUS "Defining test ${test}")
        get_filename_component(testName ${test} NAME_WE)
        add_executable(${testName} ${test})
        target_link_libraries(${testName} xolotlMemUsage)
        add_test(${testName} ${testName})
        #add a label so the tests can be run separately
        set_property(TEST ${testName} PROPERTY LABELS ${PACKAGE_NAME})
    endforeach(test ${tests})

else(Boost_FOUND)
    message(STATUS "Not building xmemusage library tests; required Boost library not available")
endif(Boost_FOUND)


