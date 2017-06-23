#ifndef XMEMUSAGE_SYSINFO_MEM_PROFILER_H
#define XMEMUSAGE_SYSINFO_MEM_PROFILER_H

#include "xolotlMemUsage/IMemUsageSampler.h"
#include "xolotlMemUsage/common/SysInfo/SysInfoSamplerBase.h"
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
std::string
GetCSVHeaderString(void)
{
    return "totalram_bytes,freeram_bytes,sharedram_bytes,bufferram_bytes,totalswap_bytes,freeswap_bytes,totalhigh_bytes,freehigh_bytes";
}

inline
std::ostream&
operator<<(std::ostream& os, const SampleAverages& sa)
{
    os << sa.totalram << ','
        << sa.freeram << ','
        << sa.sharedram << ','
        << sa.bufferram << ','
        << sa.totalswap << ','
        << sa.freeswap << ','
        << sa.totalhigh << ','
        << sa.freehigh;
    return os;
}


using ProfilerTimeHistogram = TimeHistogram<Sample, RunningSampleTotals, SampleAverages, AsyncSamplingThreadBase::ClockType>;

} // namespace SysInfo


//
// We need to forward-declare the specializations of these TimeHistogram
// methods here, or else they will be instantiated with their default
// implementations.
//
// In C++11, these explicit specializations have to occur in the
// namespace that *encloses* the specialized template.  (Earlier
// standards required them to be in the namespace "of which the template
// is a member.")  So, we close off the namespace to declare the
// specializations and later re-open it for the code that will
// cause the class to be implicitly instantiated except for these
// explicit specializations.
//

template<>
void
SysInfo::ProfilerTimeHistogram::BinData::HandleSample(const SysInfo::Sample& sample);

template<>
void
SysInfo::ProfilerTimeHistogram::BinData::Reset(void);

template<>
SysInfo::SampleAverages
SysInfo::ProfilerTimeHistogram::BinData::GetMetricValue(void) const;


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
