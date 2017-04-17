#ifndef XMEMUSAGE_STATM_MEM_PROFILER_H
#define XMEMUSAGE_STATM_MEM_PROFILER_H

#include "xolotlMemUsage/IMemUsageSampler.h"
#include "xolotlMemUsage/common/StatmSamplerBase.h"
#include "xolotlMemUsage/profile/TimeHistogram.h"

namespace xolotlMemUsage {

namespace Statm {


struct RunningSampleTotals
{
    uint64_t vmSize;
    uint64_t vmRSS;
    uint64_t rss;
    uint64_t text;
    uint64_t dataAndStack;

    RunningSampleTotals(void)
      : vmSize(0),
        vmRSS(0),
        rss(0),
        text(0),
        dataAndStack(0)
    {
        // Nothing else to do.
    }

    RunningSampleTotals& operator+=(const RunningSampleTotals& other)
    {
        vmSize += other.vmSize;
        vmRSS += other.vmRSS;
        rss += other.rss;
        text += other.text;
        dataAndStack += other.dataAndStack;

        return *this;
    }
};

struct SampleAverages
{
    double vmSize;
    double vmRSS;
    double rss;
    double text;
    double dataAndStack;

    SampleAverages(void)
      : vmSize(0),
        vmRSS(0),
        rss(0),
        text(0),
        dataAndStack(0)
    {
        // Nothing else to do.
    }
};


inline
std::ostream&
operator<<(std::ostream& os, const SampleAverages& sa)
{
    os << sa.vmSize << '\t'
        << sa.vmRSS << '\t'
        << sa.rss << '\t'
        << sa.text << '\t'
        << sa.dataAndStack;
    return os;
}


using ProfilerTimeHistogram = TimeHistogram<Sample, RunningSampleTotals, SampleAverages, AsyncSamplingThreadBase::ClockType>;


//
// For reasons I don't understand, at least some compilers (e.g., Intel 17)
// require these template specializations for structures to be done
// outside of the Statm namespace.
// We close off the namespace, define the specializations, then
// restart the namespace for the Profiler class.
//

} // namespace Statm

template<>
void
Statm::ProfilerTimeHistogram::BinData::HandleSample(const Statm::Sample& sample)
{
    nSamples++;

    runningValue.vmSize += sample.vmSize;
    runningValue.vmRSS += sample.vmRSS;
    runningValue.rss += sample.rss;
    runningValue.text += sample.text;
    runningValue.dataAndStack += sample.dataAndStack;
}

template<>
void
Statm::ProfilerTimeHistogram::BinData::Reset(void)
{
    nSamples = 0;

    runningValue.vmSize = 0;
    runningValue.vmRSS = 0;
    runningValue.rss = 0;
    runningValue.text = 0;
    runningValue.dataAndStack = 0;
}

template<>
Statm::SampleAverages
Statm::ProfilerTimeHistogram::BinData::GetMetricValue(void) const
{
    Statm::SampleAverages ret;

    if(nSamples > 0)
    {
        ret.vmSize = static_cast<double>(runningValue.vmSize) / nSamples;
        ret.vmRSS = static_cast<double>(runningValue.vmRSS) / nSamples;
        ret.rss = static_cast<double>(runningValue.rss) / nSamples;
        ret.text = static_cast<double>(runningValue.text) / nSamples;
        ret.dataAndStack = static_cast<double>(runningValue.dataAndStack) / nSamples;
    }
    return ret;
}


// Restarting Statm namespace - see comment above.
namespace Statm {

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




class Profiler : public Statm::SamplerBase<RunningProfileData>
{
private:
    static const unsigned int nTimeHistogramBins;
    static const AsyncSamplingThreadBase::ClockType::duration initialBinWidth;

public:
    Profiler(std::string _name)
      : Statm::SamplerBase<RunningProfileData>(_name,
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

} // namespace Statm

} // namespace xolotlMemUsage

#endif // XMEMUSAGE_STATM_MEM_PROFILER_H
