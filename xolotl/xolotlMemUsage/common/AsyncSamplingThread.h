#ifndef XMEMUSAGE_ASYNC_SAMPLING_THREAD_H
#define XMEMUSAGE_ASYNC_SAMPLING_THREAD_H

#include <iostream>
#include <string>
#include <tuple>
#include <future>
#include <chrono>
#include <thread>
#include <mutex>
#include <map>
#include <set>
#include <cassert>


namespace xolotlMemUsage {


class AsyncSamplingThreadBase
{
public:
    using ClockType = std::chrono::system_clock;

protected:
    static ClockType::duration samplingInterval;

public:
    static void SetSamplingInterval(ClockType::duration _interval)
    {
        samplingInterval = _interval;
    }

    AsyncSamplingThreadBase(void) = default;
};



template<typename SampleType, typename SamplingSupportType> class AsyncSampler;

template<typename SampleType, typename SamplingSupportType>
class AsyncSamplingThread : public AsyncSamplingThreadBase
{
private:

    std::mutex mtx;     // For ensuring thread-safe access to our data
    std::map<std::string, AsyncSampler<SampleType, SamplingSupportType>*> allSamplers;
    std::set<AsyncSampler<SampleType, SamplingSupportType>*> activeSamplers;
    std::future<void> samplingThreadResult;
    bool haveSamplingThread;
    bool samplingThreadDone;
    SamplingSupportType supportData;

    // Collect a sample.
    std::tuple<AsyncSamplingThreadBase::ClockType::time_point, SampleType>
        CollectSample(const SamplingSupportType& supportData) const;
    
    // Main function for asynchronous sampling thread.
    void BeSamplingThread(void);

public:
    AsyncSamplingThread(void)
      : haveSamplingThread(false),
        samplingThreadDone(true)
    {
        // Nothing else to do.
    }

    void StartSamplingFor(AsyncSampler<SampleType, SamplingSupportType>* sampler);
    void StopSamplingFor(AsyncSampler<SampleType, SamplingSupportType>* sampler);
};

} // end namespace xolotlMemUsage

#endif // XMEMUSAGE_ASYNC_SAMPLING_THREAD_H
