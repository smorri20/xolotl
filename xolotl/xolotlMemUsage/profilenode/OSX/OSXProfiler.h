#ifndef XMEMUSAGE_OSX_MEM_PROFILER_H
#define XMEMUSAGE_OSX_MEM_PROFILER_H

#include "xolotlMemUsage/IMemUsageSampler.h"
#include "xolotlMemUsage/common/OSX/OSXSamplerBase.h"
#include "xolotlMemUsage/common/TimeHistogram.h"

namespace xolotlMemUsage {

namespace OSX {

struct RunningSampleTotals
{
    double free_count;
    double active_count;
    double inactive_count;
    double wire_count;
    double zero_fill_count;
    double reactivations;
    double pageins;
    double pageouts;
    double faults;
    double cow_faults;
    double lookups;
    double hits;
    double purges;
    double purgeable_count;

    RunningSampleTotals(void)
      : free_count(0),
        active_count(0),
        inactive_count(0),
        wire_count(0),
        zero_fill_count(0),
        reactivations(0),
        pageins(0),
        pageouts(0),
        faults(0),
        cow_faults(0),
        lookups(0),
        hits(0),
        purges(0),
        purgeable_count(0)
    {
        // Nothing else to do.
    }

    RunningSampleTotals& operator+=(const RunningSampleTotals& other)
    {
        free_count += other.free_count;
        active_count += other.active_count;
        inactive_count += other.inactive_count;
        wire_count += other.wire_count;
        zero_fill_count += other.zero_fill_count;
        reactivations += other.reactivations;
        pageins += other.pageins;
        pageouts += other.pageouts;
        faults += other.faults;
        cow_faults += other.cow_faults;
        lookups += other.lookups;
        hits += other.hits;
        purges += other.purges;
        purgeable_count += other.purgeable_count;

        return *this;
    }
};

struct SampleAverages
{
    double free_count;
    double active_count;
    double inactive_count;
    double wire_count;
    double zero_fill_count;
    double reactivations;
    double pageins;
    double pageouts;
    double faults;
    double cow_faults;
    double lookups;
    double hits;
    double purges;
    double purgeable_count;


    SampleAverages(void)
      : free_count(0),
        active_count(0),
        inactive_count(0),
        wire_count(0),
        zero_fill_count(0),
        reactivations(0),
        pageins(0),
        pageouts(0),
        faults(0),
        cow_faults(0),
        lookups(0),
        hits(0),
        purges(0),
        purgeable_count(0)
    {
        // Nothing else to do.
    }
};



inline
std::string
GetCSVHeaderString(void)
{
    return "free_count,active_count,inactive_count,wire_count,zero_fill_count,reactivations,pageins,pageouts,faults,cow_faults,lookups,hits,purges,purgeable_count";
}


inline
std::ostream&
operator<<(std::ostream& os, const SampleAverages& sa)
{
    os << sa.free_count << ','
        << sa.active_count << ','
        << sa.inactive_count << ','
        << sa.wire_count << ','
        << sa.zero_fill_count << ','
        << sa.reactivations << ','
        << sa.pageins << ','
        << sa.pageouts << ','
        << sa.faults << ','
        << sa.cow_faults << ','
        << sa.lookups << ','
        << sa.hits << ','
        << sa.purges << ','
        << sa.purgeable_count;

    return os;
}


using ProfilerTimeHistogram = TimeHistogram<Sample, RunningSampleTotals, SampleAverages, AsyncSamplingThreadBase::ClockType>;


} // namespace OSX


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
OSX::ProfilerTimeHistogram::BinData::HandleSample(const OSX::Sample& sample);

template<>
void
OSX::ProfilerTimeHistogram::BinData::Reset(void);

template<>
OSX::SampleAverages
OSX::ProfilerTimeHistogram::BinData::GetMetricValue(void) const;


// Restarting OSX namespace - see comment above.
namespace OSX {

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




class Profiler : public OSX::SamplerBase<RunningProfileData>
{
private:
    static const unsigned int nTimeHistogramBins;
    static const AsyncSamplingThreadBase::ClockType::duration initialBinWidth;

public:
    Profiler(std::string _name)
      : OSX::SamplerBase<RunningProfileData>(_name,
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

} // namespace OSX

} // namespace xolotlMemUsage

#endif // XMEMUSAGE_OSX_MEM_PROFILER_H
