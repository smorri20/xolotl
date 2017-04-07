#ifndef ASYNC_SAMPLING_THREAD_H
#define ASYNC_SAMPLING_THREAD_H

#include <iostream>
#include <string>
#include <future>
#include <chrono>
#include <thread>
#include <mutex>
#include <map>
#include <set>
#include <cassert>


namespace xolotlPerf {

template<typename RPD, typename PSD> class AsyncSamplingRegion;

template<typename RPD, typename PSD>
class AsyncSamplingThread
{
public:
    typedef std::chrono::duration<double, std::milli> SamplingIntervalType;

private:
    std::mutex mtx;     // For ensuring thread-safe access to our data
    PSD supportData;    // Data needed when sampling.
    std::map<std::string, AsyncSamplingRegion<RPD, PSD>*> allRegions;
    std::set<AsyncSamplingRegion<RPD, PSD>*> activeRegions;
    std::future<void> samplingThreadResult;
    bool haveSamplingThread;
    bool samplingThreadDone;
    SamplingIntervalType samplingInterval;

    // Collect one sample, attributing it to any active profilers.
    void Sample(void);

    // Thread main program for asynchronous sampling thread.
    void CollectSamples(void);

public:
    AsyncSamplingThread(const SamplingIntervalType& _interval)
      : haveSamplingThread(false),
        samplingThreadDone(true),
        samplingInterval(_interval)
    {
        // Nothing else to do.
    }

    void StartSamplingFor(AsyncSamplingRegion<RPD, PSD>* region);
    void StopSamplingFor(AsyncSamplingRegion<RPD, PSD>* region);
};

} // end namespace xolotlPerf

#endif // ASYNC_SAMPLING_THREAD_H
