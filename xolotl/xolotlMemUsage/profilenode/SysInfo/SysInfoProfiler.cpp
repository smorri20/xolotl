#include "xolotlMemUsage/profilenode/SysInfo/SysInfoProfiler.h"
#include "xolotlMemUsage/profilenode/NodeMemUsageProfile.h"

namespace xolotlMemUsage {

namespace SysInfo {

const unsigned int Profiler::nTimeHistogramBins = 40;

const AsyncSamplingThreadBase::ClockType::duration Profiler::initialBinWidth = std::chrono::duration<uint64_t, std::milli>(1000);

std::shared_ptr<IMemUsageSampler::MemUsageData>
Profiler::GetCurrentProfile(void) const
{
    return std::make_shared<NodeMemUsageProfile>(GetRunningSampleData().hist);
}

} // end namespace SysInfo


// Define our explicit specializations.
// We do this outside the OSX namespace, as per C++11 requirement.
// (See comment in the associated header file.)
template<>
void
SysInfo::ProfilerTimeHistogram::BinData::HandleSample(const SysInfo::Sample& sample)
{
    nSamples++;

    if(runningValue.mem_unit == 0)
    {
        runningValue.mem_unit = sample.mem_unit;
    }
    else
    {
        assert(sample.mem_unit == runningValue.mem_unit);
    }

    runningValue.totalram += sample.totalram;
    runningValue.freeram += sample.freeram;
    runningValue.sharedram += sample.sharedram;
    runningValue.bufferram += sample.bufferram;
    runningValue.totalswap += sample.totalswap;
    runningValue.freeswap += sample.freeswap;
    runningValue.totalhigh += sample.totalhigh;
    runningValue.freehigh += sample.freehigh;
}

template<>
void
SysInfo::ProfilerTimeHistogram::BinData::Reset(void)
{
    nSamples = 0;

    runningValue.mem_unit = 0;

    runningValue.totalram = 0;
    runningValue.freeram = 0;
    runningValue.sharedram = 0;
    runningValue.bufferram = 0;
    runningValue.totalswap = 0;
    runningValue.freeswap = 0;
    runningValue.totalhigh = 0;
    runningValue.freehigh = 0;
}

template<>
SysInfo::SampleAverages
SysInfo::ProfilerTimeHistogram::BinData::GetMetricValue(void) const
{
    SysInfo::SampleAverages ret;

    if(nSamples > 0)
    {
        ret.totalram = runningValue.totalram / nSamples;
        ret.freeram = runningValue.freeram / nSamples;
        ret.sharedram = runningValue.sharedram / nSamples;
        ret.bufferram = runningValue.bufferram / nSamples;
        ret.totalswap = runningValue.totalswap / nSamples;
        ret.freeswap = runningValue.freeswap / nSamples;
        ret.totalhigh = runningValue.totalhigh / nSamples;
        ret.freehigh = runningValue.freehigh / nSamples;
    }
    return ret;
}

} // namespace xolotlMemUsage

