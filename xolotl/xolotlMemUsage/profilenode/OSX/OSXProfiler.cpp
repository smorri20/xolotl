#include "xolotlMemUsage/profilenode/OSX/OSXProfiler.h"
#include "xolotlMemUsage/profilenode/NodeMemUsageProfile.h"


namespace xolotlMemUsage {

namespace OSX {

const unsigned int Profiler::nTimeHistogramBins = 40;

const AsyncSamplingThreadBase::ClockType::duration Profiler::initialBinWidth = std::chrono::duration<uint64_t, std::milli>(1000);

std::shared_ptr<IMemUsageSampler::MemUsageData>
Profiler::GetCurrentProfile(void) const
{
    return std::make_shared<NodeMemUsageProfile>(GetRunningSampleData().hist);
}

} // namespace OSX


// Define our explicit specializations.
// We do this outside the OSX namespace, as per C++11 requirement.
// (See comment in the associated header file.)
template<>
void
OSX::ProfilerTimeHistogram::BinData::HandleSample(const OSX::Sample& sample)
{
    nSamples++;

    runningValue.free_count += sample.stats.free_count;
    runningValue.active_count += sample.stats.active_count;
    runningValue.inactive_count += sample.stats.inactive_count;
    runningValue.wire_count += sample.stats.wire_count;
    runningValue.zero_fill_count += sample.stats.zero_fill_count;
    runningValue.reactivations += sample.stats.reactivations;
    runningValue.pageins += sample.stats.pageins;
    runningValue.pageouts += sample.stats.pageouts;
    runningValue.faults += sample.stats.faults;
    runningValue.cow_faults += sample.stats.cow_faults;
    runningValue.lookups += sample.stats.lookups;
    runningValue.hits += sample.stats.hits;
    runningValue.purges += sample.stats.purges;
    runningValue.purgeable_count += sample.stats.purgeable_count;
}


template<>
void
OSX::ProfilerTimeHistogram::BinData::Reset(void)
{
    nSamples = 0;

    runningValue.free_count = 0;
    runningValue.active_count = 0;
    runningValue.inactive_count = 0;
    runningValue.wire_count = 0;
    runningValue.zero_fill_count = 0;
    runningValue.reactivations = 0;
    runningValue.pageins = 0;
    runningValue.pageouts = 0;
    runningValue.faults = 0;
    runningValue.cow_faults = 0;
    runningValue.lookups = 0;
    runningValue.hits = 0;
    runningValue.purges = 0;
    runningValue.purgeable_count = 0;
}

template<>
OSX::SampleAverages
OSX::ProfilerTimeHistogram::BinData::GetMetricValue(void) const
{
    OSX::SampleAverages ret;

    if(nSamples > 0)
    {
        ret.free_count = runningValue.free_count / nSamples;
        ret.active_count = runningValue.active_count / nSamples;
        ret.inactive_count = runningValue.inactive_count / nSamples;
        ret.wire_count = runningValue.wire_count / nSamples;
        ret.zero_fill_count = runningValue.zero_fill_count / nSamples;
        ret.reactivations = runningValue.reactivations / nSamples;
        ret.pageins = runningValue.pageins / nSamples;
        ret.pageouts = runningValue.pageouts / nSamples;
        ret.faults = runningValue.faults / nSamples;
        ret.cow_faults = runningValue.cow_faults / nSamples;
        ret.lookups = runningValue.lookups / nSamples;
        ret.hits = runningValue.hits / nSamples;
        ret.purges = runningValue.purges / nSamples;
        ret.purgeable_count = runningValue.purgeable_count / nSamples;
    }
    return ret;
}

} // namespace xolotlMemUsage


