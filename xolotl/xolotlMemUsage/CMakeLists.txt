#Set the package name
SET(PACKAGE_NAME "xolotl.memusage")
#Set the description
SET(PACKAGE_DESCRIPTION "Xolotl Memory Usage")
#Set the library name
SET(LIBRARY_NAME "xolotlMemUsage")

# We will be looking for specific headers.
include(CheckIncludeFiles)

#Include headers in order to build
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_BINARY_DIR})

# Specify files we always need to build/depend on, regardless of anything
# else that is available/supported.
# Note: we purposefully avoid using file(GLOB ...) here.  Since cmake
# can't tell when files are added to a directory, explicitly naming
# files follows principle of least surprise.
list(APPEND HEADERS IHandlerRegistry.h
    IMemUsageSampler.h
    RuntimeError.h
    xolotlMemUsage.h
    common/AsyncSampler.h
    common/AsyncSamplingThread.h
    common/AsyncSamplingThreadDefs.h
    common/CommonHandlerRegistry.h
    common/MemUsageSamplerBase.h
    common/TimeHistogram.h
    dummy/DummyHandlerRegistry.h
    dummy/DummyMemUsageSampler.h)
list(APPEND SRCS xolotlMemUsage.cpp
    common/AsyncSamplingThreadBase.cpp
    common/CommonHandlerRegistry.cpp
    common/TimeHistogram.cpp
    dummy/DummyHandlerRegistry.cpp)

# Check if we have sys/sysinfo.h.  We use this on Linux for per-node data.
check_include_files(sys/sysinfo.h HAVE_SYSINFO)
if(HAVE_SYSINFO)
    set(HAVE_PER_NODE_DATA_SOURCE 1)
    list(APPEND HEADERS common/SysInfo/SysInfoSamplerBase.h
        common/SysInfo/SysInfoSamplingThread.h
        summarynode/SysInfo/SysInfoSampler.h
        summarynode/SysInfo/NodeMemUsageStats.h
        profilenode/SysInfo/SysInfoProfiler.h)
    list(APPEND SRCS common/SysInfo/SysInfoSamplingThread.cpp
        summarynode/SysInfo/SysInfoSampler.cpp
        summarynode/SysInfo/NodeMemUsageStats.cpp
        profilenode/SysInfo/SysInfoProfiler.cpp)
endif()

# Check if we have support for node-level host_statistics (i.e., OS X).
check_include_files("mach/host_info.h;mach/mach_host.h;mach/task_info.h;mach/task.h" HAVE_MACH_HOST_STATISTICS)
if(HAVE_MACH_HOST_STATISTICS)
    set(HAVE_PER_NODE_DATA_SOURCE 1)
    list(APPEND HEADERS common/OSX/OSXSamplerBase.h
        common/OSX/OSXSamplingThread.h
        summarynode/OSX/OSXSampler.h
        summarynode/OSX/NodeMemUsageStats.h
        profilenode/OSX/OSXInfoProfiler.h)
    list(APPEND SRCS common/OSX/OSXSamplingThread.cpp
        summarynode/OSX/OSXSampler.cpp
        summarynode/OSX/NodeMemUsageStats.cpp
        profilenode/OSX/OSXProfiler.cpp)
endif()

if(HAVE_PER_NODE_DATA_SOURCE)
    list(APPEND HEADERS summarynode/NodeMemUsageObjStatistics.h
        summarynode/NodeMemUsageSampler.h
        summarynode/NodeMemUsageStatsBase.h
        summarynode/SummaryNodeHandlerRegistry.h
        profilenode/NodeMemUsageProfile.h
        profilenode/NodeMemUsageProfiler.h
        profilenode/ProfileNodeHandlerRegistry.h)
    list(APPEND SRCS summarynode/SummaryNodeHandlerRegistry.cpp
        profilenode/MemUsageProfile.cpp
        profilenode/ProfileNodeHandlerRegistry.cpp)
endif()

# Check if we have /proc.  We assume that if we do, we will have
# /proc/<pid>/statm for our per-process data source.
message(STATUS "Looking for /proc")
if(EXISTS "/proc")
    message(STATUS "Looking for /proc - found")
    set(HAVE_STATM TRUE)
else()
    message(STATUS "Looking for /proc - not found")
    set(HAVE_STATM FALSE)
endif()

if(HAVE_STATM)
    set(HAVE_PER_PROC_DATA_SOURCE 1)
    list(APPEND HEADERS common/Statm/StatmSamplerBase.h
        common/Statm/StatmSamplingThread.h
        summaryproc/Statm/StatmSampler.h
        summaryproc/Statm/MemUsageStats.h
        profileproc/Statm/StatmProfiler.h)
    list(APPEND SRCS common/Statm/StatmSamplingThread.cpp
        summaryproc/Statm/StatmSampler.cpp
        summaryproc/Statm/MemUsageStats.cpp
        profileproc/Statm/StatmProfiler.cpp)
endif()

if(HAVE_PER_PROC_DATA_SOURCE)
    list(APPEND HEADERS summaryproc/MemUsageObjStatistics.h
        summaryproc/MemUsageSampler.h
        summaryproc/MemUsageStatsBase.h
        summaryproc/SummaryProcHandlerRegistry.h
        profileproc/MemUsageProfile.h
        profileproc/MemUsageProfiler.h
        profileproc/ProfileProcHandlerRegistry.h)
    list(APPEND SRCS summaryproc/SummaryProcHandlerRegistry.cpp
        profileproc/MemUsageProfile.cpp
        profileproc/ProfileProcHandlerRegistry.cpp)
endif()

# Specify the library to build
add_library(${LIBRARY_NAME} STATIC ${SRCS})

#Install the xolotl header files
install(FILES ${HEADERS} DESTINATION include)
install(TARGETS ${LIBRARY_NAME} DESTINATION lib)

# Configure file for the memory usage library.
# Note: must do this after all the other checks, or else the 
# contents of the generated file won't take the check results into account.
configure_file( "${CMAKE_CURRENT_SOURCE_DIR}/memUsageConfig.h.in" "${CMAKE_CURRENT_BINARY_DIR}/memUsageConfig.h")


