#ifndef XMEMUSAGE_NODE_MEM_USAGE_NODE_PROFILER_H
#define XMEMUSAGE_NODE_MEM_USAGE_NODE_PROFILER_H

#include "xolotlMemUsage/common/MemUsageSamplerBase.h"
#include "xolotlMemUsage/profilenode/SysInfoProfiler.h"


namespace xolotlMemUsage {

class NodeMemUsageProfiler : public MemUsageSamplerBase<SysInfo::Profiler> {

public:

    /// Disallow building an object without a name.
    NodeMemUsageProfiler(void) = delete;


    /**
     * Construct a MemUsageProfiler with given name.
     *
     * @param name The object's name.
     */
    NodeMemUsageProfiler(const std::string& name)
      : MemUsageSamplerBase<SysInfo::Profiler>(name) {

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

#endif // XMEMUSAGE_NODE_MEM_USAGE_NODE_PROFILER_H
