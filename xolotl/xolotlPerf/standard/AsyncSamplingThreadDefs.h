#ifndef ASYNC_SAMPLING_THREAD_DEFS_H
#define ASYNC_SAMPLING_THREAD_DEFS_H

// This file contains implementations of AsyncSamplingThread template 
// class member functions.  It is needed to break a dependency cycle
// between AsyncSamplingRegion and AsyncSamplingThread.
#include "AsyncSamplingThread.h"

namespace xolotlPerf {

template<typename RPD, typename PSD>
void
AsyncSamplingThread<RPD, PSD>::StartSamplingFor(AsyncSamplingRegion<RPD, PSD>* region)
{
    // Ensure we are the only ones updating our data structures.
    std::lock_guard<std::mutex> ourLock(mtx);

    // Ensure we have a sampling thread.
    if(not haveSamplingThread)
    {
        // We do not have a sampling thread.  Create one.
        haveSamplingThread = true;
        samplingThreadDone = false;
        samplingThreadResult = std::async(&AsyncSamplingThread::CollectSamples, this);
    }

    // Ensure we know about the region.
    auto arIter = allRegions.find(region->GetName());
    if(arIter == allRegions.end())
    {
        // We did not already know about a region with this name.
        allRegions.emplace(region->GetName(), region);
    }
    else if(arIter->second != region)
    {
        // We knew about a region with this name, but it isn't 
        // the one that we have been given.
        std::cerr << "WARNING: unmatched start/stop sampling for async sampling regions."
            << "Replacing region with name " << region->GetName() 
            << ".  Old region's sampled data will be lost." 
            << std::endl;
        auto oldRegion = arIter->second;
        allRegions[region->GetName()] = region;

        // Ensure that the old region is no longer considered active.
        auto activeIter = activeRegions.find(oldRegion);
        if(activeIter != activeRegions.end())
        {
            // Old region object was active, so remove it so that new
            // region object can safely take its place.
            activeRegions.erase(activeIter);
        }
    }

    // Ensure the region is considered active.
    activeRegions.emplace(region);
}


template<typename RPD, typename PSD>
void
AsyncSamplingThread<RPD, PSD>::StopSamplingFor(AsyncSamplingRegion<RPD, PSD>* region)
{
    // Ensure we are the only ones updating our data structures.
    std::lock_guard<std::mutex> ourLock(mtx);

    // Remove the region from the set of active sampling regions.
    auto activeIter = activeRegions.find(region);
    if(activeIter != activeRegions.end())
    {
        activeRegions.erase(activeIter);

        // Check whether we have *any* regions that are active.
        if(activeRegions.empty())
        {
            // We have no active sampling regions - no need
            // to keep the sampling thread running.
            samplingThreadDone = true;        
            samplingThreadResult.get();   // wait till sample thread goes away
            haveSamplingThread = false;
        }
    }
}

template<typename RPD, typename PSD>
void
AsyncSamplingThread<RPD, PSD>::CollectSamples(void)
{
    while(not samplingThreadDone)
    {
        Sample();
        std::this_thread::sleep_for(samplingInterval);
    }
}

} // end namespace xolotlPerf

#endif // ASYNC_SAMPLING_THREAD_DEFS_H
