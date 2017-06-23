#ifndef XMEMUSAGE_OSX_MEM_SAMPLER_H
#define XMEMUSAGE_OSX_MEM_SAMPLER_H

#include "xolotlMemUsage/IMemUsageSampler.h"
#include "xolotlMemUsage/common/OSX/OSXSamplerBase.h"

namespace xolotlMemUsage {

namespace OSX {

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

    RunningMetricStats free_count;
    RunningMetricStats active_count;
    RunningMetricStats inactive_count;
    RunningMetricStats wire_count;
    RunningMetricStats zero_fill_count;
    RunningMetricStats reactivations;
    RunningMetricStats pageins;
    RunningMetricStats pageouts;
    RunningMetricStats faults;
    RunningMetricStats cow_faults;
    RunningMetricStats lookups;
    RunningMetricStats hits;
    RunningMetricStats purges;
    RunningMetricStats purgeable_count;

    RunningData(void)
      : nSamples(0)
    {
        // Nothing else to do.
    }

    void HandleNewSample(const Sample& sample)
    {
        free_count.HandleNewSample(sample.stats.free_count);
        active_count.HandleNewSample(sample.stats.active_count);
        inactive_count.HandleNewSample(sample.stats.inactive_count);
        wire_count.HandleNewSample(sample.stats.wire_count);
        zero_fill_count.HandleNewSample(sample.stats.zero_fill_count);
        reactivations.HandleNewSample(sample.stats.reactivations);
        pageins.HandleNewSample(sample.stats.pageins);
        pageouts.HandleNewSample(sample.stats.pageouts);
        faults.HandleNewSample(sample.stats.faults);
        cow_faults.HandleNewSample(sample.stats.cow_faults);
        lookups.HandleNewSample(sample.stats.lookups);
        hits.HandleNewSample(sample.stats.hits);
        purges.HandleNewSample(sample.stats.purges);
        purgeable_count.HandleNewSample(sample.stats.purgeable_count);

        nSamples++;
    }
};


class Sampler : public OSX::SamplerBase<RunningData>
{
public:
    Sampler(std::string _name)
      : OSX::SamplerBase<RunningData>(_name)
    {
        // Nothing else to do.
    }

    virtual void HandleNewSample(AsyncSamplingThreadBase::ClockType::time_point /* timestamp */, const Sample& sample)
    {
        runningData.HandleNewSample(sample);
    }

    std::shared_ptr<IMemUsageSampler::MemUsageData> GetCurrentStats(void) const;
};

} // end namespace OSX

} // end namespace xolotlMemUsage

#endif // XMEMUSAGE_OSX_MEM_SAMPLER_H
