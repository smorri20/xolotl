#ifndef XMEMUSAGE_MEM_USAGE_PROFILER_H
#define XMEMUSAGE_MEM_USAGE_PROFILER_H

#include "xolotlMemUsage/common/MemUsageSamplerBase.h"
#include "xolotlMemUsage/profile/StatmProfiler.h"


namespace xolotlMemUsage {

class MemUsageProfiler : public MemUsageSamplerBase<Statm::Profiler> {

public:

    /// Disallow building an object without a name.
    MemUsageProfiler(void) = delete;


    /**
     * Construct a MemUsageProfiler with given name.
     *
     * @param name The object's name.
     */
    MemUsageProfiler(const std::string& name)
      : MemUsageSamplerBase<Statm::Profiler>(name) {

        // Nothing else to do.
    }

    /**
     * Obtain the data we've collected.
     */
    virtual std::shared_ptr<IMemUsageSampler::MemUsageData> getValue(void) const {
        return GetCurrentProfile();
    }
};

} // namespace xolotlMemUsage

#endif // XMEMUSAGE_MEM_USAGE_PROFILER_H
