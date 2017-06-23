#ifndef XMEMUSAGE_NODE_MEM_USAGE_NODE_PROFILER_H
#define XMEMUSAGE_NODE_MEM_USAGE_NODE_PROFILER_H

#include "xolotlMemUsage/memUsageConfig.h"
#include "xolotlMemUsage/common/MemUsageSamplerBase.h"

#if defined(HAVE_SYSINFO)
    #include "xolotlMemUsage/profilenode/SysInfo/SysInfoProfiler.h"
    namespace PerNodeDataSource = xolotlMemUsage::SysInfo;
#elif defined(HAVE_MACH_HOST_STATISTICS)
    #include "xolotlMemUsage/profilenode/OSX/OSXProfiler.h"
    namespace PerNodeDataSource = xolotlMemUsage::OSX;
#else
    #error "Configuration error: thought we had a per-node data source, but no actual data source available."
#endif // defined(HAVE_SYSINFO)



namespace xolotlMemUsage {

class NodeMemUsageProfiler : public MemUsageSamplerBase<PerNodeDataSource::Profiler> {

public:

    /// Disallow building an object without a name.
    NodeMemUsageProfiler(void) = delete;


    /**
     * Construct a MemUsageProfiler with given name.
     *
     * @param name The object's name.
     */
    NodeMemUsageProfiler(const std::string& name)
      : MemUsageSamplerBase<PerNodeDataSource::Profiler>(name) {

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
