#ifndef XMEMUSAGE_SYSINFO_MEM_PROFILER_H
#define XMEMUSAGE_SYSINFO_MEM_PROFILER_H

#include "xolotlMemUsage/IMemUsageSampler.h"
#include "xolotlMemUsage/common/SysInfoSamplerBase.h"
#include "xolotlMemUsage/common/TimeHistogram.h"

namespace xolotlMemUsage {

namespace SysInfo {


struct RunningSampleTotals
{
    double totalram;
    double freeram;
    double sharedram;
    double bufferram;
    double totalswap;
    double freeswap;
    double totalhigh;
    double freehigh;
    uint32_t mem_unit;

    RunningSampleTotals(void)
      : totalram(0),
        freeram(0),
        sharedram(0),
        bufferram(0),
        totalswap(0),
        freeswap(0),
        totalhigh(0),
        freehigh(0),
        mem_unit(0)
    {
        // Nothing else to do.
    }

    RunningSampleTotals& operator+=(const RunningSampleTotals& other)
    {
        if(mem_unit == 0)
        {
            mem_unit = other.mem_unit;
        }
        else
        {
            assert(other.mem_unit == mem_unit);
        }

        totalram += other.totalram;
        freeram += other.freeram;
        sharedram += other.sharedram;
        bufferram += other.bufferram;
        totalswap += other.totalswap;
        freeswap += other.freeswap;
        totalhigh += other.totalhigh;
        freehigh += other.freehigh;

        return *this;
    }
};

struct SampleAverages
{
    double totalram;
    double freeram;
    double sharedram;
    double bufferram;
    double totalswap;
    double freeswap;
    double totalhigh;
    double freehigh;
    uint32_t mem_unit;


    SampleAverages(void)
      : totalram(0),
        freeram(0),
        sharedram(0),
        bufferram(0),
        totalswap(0),
        freeswap(0),
        totalhigh(0),
        freehigh(0),
        mem_unit(0)
    {
        // Nothing else to do.
    }
};


inline
std::ostream&
operator<<(std::ostream& os, const SampleAverages& sa)
{
    os << sa.totalram << '\t'
        << sa.freeram << '\t'
        << sa.sharedram << '\t'
        << sa.bufferram << '\t'
        << sa.totalswap << '\t'
        << sa.freeswap << '\t'
        << sa.totalhigh << '\t'
        << sa.freehigh;
    return os;
}


using ProfilerTimeHistogram = TimeHistogram<Sample, RunningSampleTotals, SampleAverages, AsyncSamplingThreadBase::ClockType>;


//
// For reasons I don't understand, at least some compilers (e.g., Intel 17)
// require these template specializations for structures to be done
// outside of the SysInfo namespace.
// We close off the namespace, define the specializations, then
// restart the namespace for the Profiler class.
//

} // namespace SysInfo

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


// Restarting SysInfo namespace - see comment above.
namespace SysInfo {

struct RunningProfileData
{
    ProfilerTimeHistogram hist;


    /// Disallow construction without arguments for constructing the histogram.
    RunningProfileData() = delete;


    RunningProfileData(unsigned int _nBins,
                        AsyncSamplingThreadBase::ClockType::duration _initialBinWidth,
                        AsyncSamplingThreadBase::ClockType::time_point _initialTimestamp)
      : hist(_nBins, _initialBinWidth, _initialTimestamp)
    {
        // Nothing else to do.
    }


    void HandleNewSample(AsyncSamplingThreadBase::ClockType::time_point timestamp, const Sample& sample)
    {
        hist.HandleSample(timestamp, sample);
    }
};




class Profiler : public SysInfo::SamplerBase<RunningProfileData>
{
private:
    static const unsigned int nTimeHistogramBins;
    static const AsyncSamplingThreadBase::ClockType::duration initialBinWidth;

public:
    Profiler(std::string _name)
      : SysInfo::SamplerBase<RunningProfileData>(_name,
                    RunningProfileData(nTimeHistogramBins, initialBinWidth, AsyncSamplingThreadBase::ClockType::now()))
    {
        // Nothing else to do.
    }

    virtual void HandleNewSample(AsyncSamplingThreadBase::ClockType::time_point timestamp, const Sample& sample)
    {
        runningData.HandleNewSample(timestamp, sample);
    }

    std::shared_ptr<IMemUsageSampler::MemUsageData> GetCurrentProfile(void) const;
};

} // namespace SysInfo

} // namespace xolotlMemUsage

#endif // XMEMUSAGE_SYSINFO_MEM_PROFILER_H
