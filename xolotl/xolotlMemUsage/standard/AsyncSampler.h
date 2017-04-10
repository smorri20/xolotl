#ifndef XMEMUSAGE_ASYNC_SAMPLER_H
#define XMEMUSAGE_ASYNC_SAMPLER_H

#include "AsyncSamplingThread.h"

namespace xolotlMemUsage {


template<typename RPD, typename PSD>
class AsyncSampler
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
    AsyncSampler(std::string _name)
      : name(_name)
    {
        // Nothing else to do.
    }

    ~AsyncSampler(void)
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

} // namespace xolotlMemUsage

#include "AsyncSamplingThreadDefs.h"

#endif // XMEMUSAGE_ASYNC_SAMPLER_H
