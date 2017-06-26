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
    list(APPEND tests DummyMemUsageSamplerTester.cpp
                        DummyHandlerRegistryTester.cpp)

    # Check if we have a usable per-proc data source.
    include(CheckMemUsageDataSources)    
    check_memusage_per_proc_data_source(HAVE_PER_PROC_DATA_SOURCE HAVE_STATM)
    if(HAVE_PER_PROC_DATA_SOURCE)
#        list(APPEND tests MemUsageSamplerTester.cpp
#                            SummaryProcHandlerRegistryTester.cpp
#                            MemUsageProfilerTester.cpp
#                            ProfileProcHandlerRegistryTester.cpp
#                            )
    endif()

    # Check if we have a usable per-node data source.
    check_memusage_per_node_data_source(HAVE_PER_NODE_DATA_SOURCE HAVE_SYSINFO HAVE_MACH_HOST_STATISTICS)
    if(HAVE_PER_NODE_DATA_SOURCE)
        list(APPEND tests SummaryNodeHandlerRegistryTester.cpp
                            ProfileNodeHandlerRegistryTester.cpp)
    endif()

    # Our profile tests can use a filesystem library (C++ std, or Boost).
    find_package(Boost COMPONENTS filesystem)
    if(Boost_FILESYSTEM_FOUND)
        set(HAVE_BOOST_FILESYSTEM TRUE)
    endif()

    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/memUsageTestConfig.h.in" "${CMAKE_CURRENT_BINARY_DIR}/memUsageTestConfig.h")

    # Define executables for each of the testers.
    foreach(test ${tests})
        message(STATUS "Defining test ${test}")
        get_filename_component(testName ${test} NAME_WE)
        add_executable(${testName} ${test})
        target_link_libraries(${testName} xolotlMemUsage ${Boost_LIBRARIES})
        add_test(${testName} ${testName})
        #add a label so the tests can be run separately
        set_property(TEST ${testName} PROPERTY LABELS ${PACKAGE_NAME})
    endforeach(test ${tests})

else(Boost_FOUND)
    message(STATUS "Not building xmemusage library tests; required Boost library not available")
endif(Boost_FOUND)


