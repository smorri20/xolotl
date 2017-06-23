#ifndef XMEMUSAGE_NODE_MEM_USAGE_SAMPLER_H
#define XMEMUSAGE_NODE_MEM_USAGE_SAMPLER_H

#include "xolotlMemUsage/memUsageConfig.h"
#include "xolotlMemUsage/common/MemUsageSamplerBase.h"

#if defined(HAVE_SYSINFO)
    #include "xolotlMemUsage/summarynode/SysInfo/SysInfoSampler.h"
    namespace PerNodeDataSource = xolotlMemUsage::SysInfo;
#elif defined(HAVE_MACH_HOST_STATISTICS)
    #include "xolotlMemUsage/summarynode/OSX/OSXSampler.h"
    namespace PerNodeDataSource = xolotlMemUsage::OSX;
#else
    #error "Configuration error: thought we had a per-node data source, but no actual data source available."
#endif // defined(HAVE_SYSINFO)


namespace xolotlMemUsage {

class NodeMemUsageSampler : public MemUsageSamplerBase<PerNodeDataSource::Sampler> {

public:
    /**
     * Disallow building an object without a name.
     */
    NodeMemUsageSampler(void) = delete;


	/**
     * Construct a MemUsageSampler with the given name.
	 *
	 * @param name The object's name.
	 */
	NodeMemUsageSampler(const std::string& name)
      : MemUsageSamplerBase<PerNodeDataSource::Sampler>(name) {

        // Nothing else to do.
    }

	/**
	 * Obtain the sampler's metric values.
	 */
    virtual std::shared_ptr<IMemUsageSampler::MemUsageData> getValue(void) const {

        return GetCurrentStats();
    }
};

} // end namespace xolotlMemUsage

#endif // XMEMUSAGE_NODE_MEM_USAGE_SAMPLER_H
