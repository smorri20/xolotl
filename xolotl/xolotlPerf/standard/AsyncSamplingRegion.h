#ifndef ASYNC_SAMPLING_REGION_H
#define ASYNC_SAMPLING_REGION_H

#include "AsyncSamplingThread.h"

namespace xolotlPerf {


template<typename RPD, typename PSD>
class AsyncSamplingRegion
{
public:
    typedef AsyncSamplingThread<RPD, PSD> SamplingThreadType;

private:
    static typename SamplingThreadType::SamplingIntervalType samplingInterval;

    static std::shared_ptr<SamplingThreadType> GetSamplingThread(void)
    {
        // This is thread-safe in C++11, but *NOT* under earlier C++ standards.
        static std::shared_ptr<SamplingThreadType> theSamplingThread = 
            std::make_shared<SamplingThreadType>(samplingInterval);
        return theSamplingThread;
    }


    std::string name;
    RPD runningSampleData;

public:
    AsyncSamplingRegion(std::string _name)
      : name(_name)
    {
        // Nothing else to do.
    }

    ~AsyncSamplingRegion(void)
    {
        StopSampling();
    }

    void
    StartSampling(void)
    {
        GetSamplingThread()->StartSamplingFor(this);
    }

    void
    StopSampling(void)
    {
        GetSamplingThread()->StopSamplingFor(this);
    }

    static void SetSamplingInterval(const typename SamplingThreadType::SamplingIntervalType& interval)
    {
        samplingInterval = interval;
    }

    std::string GetName(void) const { return name; }

    const RPD& GetRunningSampleData(void) const { return runningSampleData; }
    RPD& GetRunningSampleData(void) { return runningSampleData; }
};

} // namespace xolotlPerf

#include "AsyncSamplingThreadDefs.h"

#endif // ASYNC_SAMPLING_REGION_H
