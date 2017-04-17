#ifndef XMEMUSAGE_ASYNC_SAMPLER_H
#define XMEMUSAGE_ASYNC_SAMPLER_H

#include "AsyncSamplingThread.h"

namespace xolotlMemUsage {


template<typename SampleType, typename SamplingSupportType>
class AsyncSampler
{
public:
    using SamplingThreadType = AsyncSamplingThread<SampleType, SamplingSupportType>;

private:
    static std::shared_ptr<SamplingThreadType> GetSamplingThread(void)
    {
        // This is thread-safe in C++11, but *NOT* under earlier C++ standards.
        static std::shared_ptr<SamplingThreadType> theSamplingThread = 
            std::make_shared<SamplingThreadType>();
        return theSamplingThread;
    }

    std::string name;

public:
    AsyncSampler(std::string _name)
      : name(_name)
    {
        // Nothing else to do.
    }

    virtual ~AsyncSampler(void)
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

    std::string GetName(void) const { return name; }

    virtual void HandleNewSample(AsyncSamplingThreadBase::ClockType::time_point timestamp, const SampleType& sample) = 0;
};

} // namespace xolotlMemUsage

#include "AsyncSamplingThreadDefs.h"

#endif // XMEMUSAGE_ASYNC_SAMPLER_H
