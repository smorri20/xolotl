

macro(CHECK_MEMUSAGE_PER_PROC_DATA_SOURCE HAVE_PER_PROC_DATA_SOURCE HAVE_STATM)

    include(CheckIncludeFiles)

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
        set(HAVE_PER_PROC_DATA_SOURCE TRUE)
    endif()

endmacro()


macro(CHECK_MEMUSAGE_PER_NODE_DATA_SOURCE HAVE_PER_NODE_DATA_SOURCE HAVE_SYSINFO HAVE_MACH_HOST_STATISTICS)

    include(CheckIncludeFiles)

    # Check if we have sys/sysinfo.h.  We use this on Linux for per-node data.
    check_include_files(sys/sysinfo.h HAVE_SYSINFO)
    if(HAVE_SYSINFO)
        set(HAVE_PER_NODE_DATA_SOURCE TRUE)
    endif()

    # Check if we have support for node-level host_statistics (i.e., OS X).
    check_include_files("mach/host_info.h;mach/mach_host.h;mach/task_info.h;mach/task.h" HAVE_MACH_HOST_STATISTICS)
    if(HAVE_MACH_HOST_STATISTICS)
        set(HAVE_PER_NODE_DATA_SOURCE TRUE)
    endif()

endmacro()


