#ifndef XMEMUSAGE_SYSINFO_MEM_SAMPLER_H
#define XMEMUSAGE_SYSINFO_MEM_SAMPLER_H

#include "xolotlMemUsage/IMemUsageSampler.h"
#include "xolotlMemUsage/common/SysInfoSamplerBase.h"

namespace xolotlMemUsage {

namespace SysInfo {

// TODO share this higher up to avoid reproducing RunningMetricStats
struct RunningMetricStats
{
    uint64_t total;
    uint64_t min;
    uint64_t max;

    RunningMetricStats(void)
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


struct MetricStats
{
    uint64_t min;
    uint64_t max;
    double average;

    MetricStats(uint64_t _min = std::numeric_limits<uint64_t>::max(),
                uint64_t _max = std::numeric_limits<uint64_t>::min(),
                double _average = std::numeric_limits<double>::quiet_NaN())
      : min(_min),
        max(_max),
        average(_average)
    {
        // Nothing else to do.
    }

    MetricStats(const RunningMetricStats& rps, uint64_t nSamples)
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

    RunningMetricStats totalram;
    RunningMetricStats freeram;
    RunningMetricStats sharedram;
    RunningMetricStats bufferram;
    RunningMetricStats totalswap;
    RunningMetricStats freeswap;
    RunningMetricStats totalhigh;
    RunningMetricStats freehigh;
    RunningMetricStats mem_unit;

    RunningData(void)
      : nSamples(0)
    {
        // Nothing else to do.
    }

    void HandleNewSample(const Sample& sample)
    {
        totalram.HandleNewSample(sample.totalram);
        freeram.HandleNewSample(sample.freeram);
        sharedram.HandleNewSample(sample.sharedram);
        bufferram.HandleNewSample(sample.bufferram);
        totalswap.HandleNewSample(sample.totalswap);
        freeswap.HandleNewSample(sample.freeswap);
        totalhigh.HandleNewSample(sample.totalhigh);
        freehigh.HandleNewSample(sample.freehigh);
        mem_unit.HandleNewSample(sample.mem_unit);

        nSamples++;
    }
};


class Sampler : public SysInfo::SamplerBase<RunningData>
{
public:
    Sampler(std::string _name)
      : SysInfo::SamplerBase<RunningData>(_name)
    {
        // Nothing else to do.
    }

    virtual void HandleNewSample(AsyncSamplingThreadBase::ClockType::time_point /* timestamp */, const Sample& sample)
    {
        runningData.HandleNewSample(sample);
    }

    std::shared_ptr<IMemUsageSampler::MemUsageData> GetCurrentStats(void) const;
};

} // end namespace SysInfo

} // end namespace xolotlMemUsage

#endif // XMEMUSAGE_SYSINFO_MEM_SAMPLER_H
