#ifndef XSTATMUSAGE_STATM_SAMPLER_BASE_H
#define XSTATMUSAGE_STATM_SAMPLER_BASE_H

#include <sstream>
#include <unistd.h>
#include "xolotlMemUsage/common/AsyncSampler.h"

namespace xolotlMemUsage {

namespace Statm {

struct Sample
{
    uint64_t vmSize;
    uint64_t vmRSS;
    uint64_t rss;
    uint64_t text;
    uint64_t dataAndStack;


    Sample(void)
      : vmSize(0),
        vmRSS(0),
        rss(0),
        text(0),
        dataAndStack(0)
    {
        // Nothing else to do.
    }
};

struct SupportData
{
    std::string statmFilePath;

    SupportData(void)
      : statmFilePath([]{std::ostringstream ostr; ostr << "/proc/" << getpid() << "/statm"; return ostr.str();}())
    {
        // Nothing else to do.
    }
};


template<typename RunningDataType>
class SamplerBase : public AsyncSampler<Sample, SupportData>
{
protected:
    RunningDataType runningData;

public:
    /// Disallow creating a sampler without a name.
    SamplerBase(void) = delete;

    /**
     * Construct a SamplerBase with the given name.
     *
     * @param _name The name to associate with the newly-created SamplerBase.
     */
    SamplerBase(std::string _name)
      : AsyncSampler(_name)
    {
        // Nothing else to do.
    }

    /**
     * Construct a SamplerBase with given name and given
     * RunningDataType as a starting point.
     *
     * @ param _name Name to associate with newly-created SamplerBase.
     * @ param _runningData Initial state of running data to use.
     */
    SamplerBase(std::string _name, const RunningDataType& _runningData)
      : AsyncSampler(_name),
        runningData(_runningData)
    {
        // Nothing else to do.
    }

    const RunningDataType& GetRunningSampleData(void) const
    {
        return runningData;
    }
};

} // namespace Statm

} // end namespace xolotlMemUsage

#endif // XSTATMUSAGE_STATM_SAMPLER_BASE_H
