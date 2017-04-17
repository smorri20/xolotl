#ifndef XSYSINFOUSAGE_SYSINFO_SAMPLER_BASE_H
#define XSYSINFOUSAGE_SYSINFO_SAMPLER_BASE_H

#include <sstream>
#include <sys/sysinfo.h>
#include "xolotlMemUsage/common/AsyncSampler.h"

namespace xolotlMemUsage {

namespace SysInfo {

struct Sample
{
    uint64_t totalram;
    uint64_t freeram;
    uint64_t sharedram;
    uint64_t bufferram;
    uint64_t totalswap;
    uint64_t freeswap;
    uint64_t totalhigh;
    uint64_t freehigh;
    uint64_t mem_unit;

    Sample(void)
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

    Sample(const struct sysinfo& si)
      : totalram(si.totalram),
        freeram(si.freeram),
        sharedram(si.sharedram),
        bufferram(si.bufferram),
        totalswap(si.totalswap),
        freeswap(si.freeswap),
        totalhigh(si.totalhigh),
        freehigh(si.freehigh),
        mem_unit(si.mem_unit)
    {
        // Nothing else to do.
    }
};


struct SupportData
{
    SupportData(void)
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

} // namespace SysInfo

} // end namespace xolotlMemUsage

#endif // XSYSINFOUSAGE_SYSINFO_SAMPLER_BASE_H
