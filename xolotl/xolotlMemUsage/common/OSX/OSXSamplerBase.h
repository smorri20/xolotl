#ifndef XOSXUSAGE_OSX_SAMPLER_BASE_H
#define XOSXUSAGE_OSX_SAMPLER_BASE_H

#include <sstream>
#include <cstring>
#include <mach/host_info.h>
#include "xolotlMemUsage/common/AsyncSampler.h"

namespace xolotlMemUsage {

namespace OSX {

struct Sample
{
    vm_statistics64_data_t stats;
    

    Sample(void)
    {
        std::memset(&stats, 0, sizeof(stats));
    }

    Sample(const vm_statistics64_data_t& _other)
      : stats(_other)
    {
        // Nothing else to do.
    }
};


struct SupportData
{
    SupportData(void)
    {
        // Nothing else to do for this type.
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

} // namespace OSX

} // end namespace xolotlMemUsage

#endif // XOSXUSAGE_OSX_SAMPLER_BASE_H
