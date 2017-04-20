#ifndef XMEMUSAGE_MEM_SAMPLER_H
#define XMEMUSAGE_MEM_SAMPLER_H

#include "xolotlMemUsage/IMemUsageSampler.h"
#include "xolotlMemUsage/common/StatmSamplerBase.h"

namespace xolotlMemUsage {

namespace Statm {

struct RunningPageStats
{
    uint64_t total;
    uint64_t min;
    uint64_t max;

    RunningPageStats(void)
      : total(0),
        min(std::numeric_limits<uint64_t>::max()),
        max(std::numeric_limits<uint64_t>::min())
    {
        // Nothing else to do.
    }

    void HandleNewSample(uint64_t val)
    {
        total += val;
        if(val < min)
        {
            min = val;
        }
        if(val > max)
        {
            max = val;
        }
    }
};


struct PageStats
{
    uint64_t min;
    uint64_t max;
    double average;

    PageStats(uint64_t _min = std::numeric_limits<uint64_t>::max(),
                uint64_t _max = std::numeric_limits<uint64_t>::min(),
                double _average = std::numeric_limits<double>::quiet_NaN())
      : min(_min),
        max(_max),
        average(_average)
    {
        // Nothing else to do.
    }

    PageStats(const RunningPageStats& rps, uint64_t nSamples)
      : min(rps.min),
        max(rps.max),
        average(((double)rps.total) / nSamples)
    {
        // Nothing else to do.
    }

};


struct RunningData
{
    uint64_t nSamples;

    RunningPageStats vmSize;
    RunningPageStats vmRSS;
    RunningPageStats rss;
    RunningPageStats text;
    RunningPageStats dataAndStack;

    RunningData(void)
      : nSamples(0)
    {
        // Nothing else to do.
    }

    void HandleNewSample(const Sample& sample)
    {
        vmSize.HandleNewSample(sample.vmSize);
        vmRSS.HandleNewSample(sample.vmRSS);
        rss.HandleNewSample(sample.rss);
        text.HandleNewSample(sample.text);
        dataAndStack.HandleNewSample(sample.dataAndStack);
        nSamples++;
    }
};


class Sampler : public Statm::SamplerBase<RunningData>
{
public:
    Sampler(std::string _name)
      : Statm::SamplerBase<RunningData>(_name)
    {
        // Nothing else to do.
    }

    virtual void HandleNewSample(AsyncSamplingThreadBase::ClockType::time_point /* timestamp */, const Sample& sample)
    {
        runningData.HandleNewSample(sample);
    }

    std::shared_ptr<IMemUsageSampler::MemUsageData> GetCurrentStats(void) const;
};

} // end namespace Statm

} // end namespace xolotlMemUsage

#endif // XMEMUSAGE_MEM_SAMPLER_H
