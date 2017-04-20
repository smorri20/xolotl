#ifndef XMEMUSAGE_NODE_MEM_USAGE_SAMPLER_H
#define XMEMUSAGE_NODE_MEM_USAGE_SAMPLER_H

#include "xolotlMemUsage/common/MemUsageSamplerBase.h"
#include "xolotlMemUsage/summarynode/SysInfoSampler.h"


namespace xolotlMemUsage {

class NodeMemUsageSampler : public MemUsageSamplerBase<SysInfo::Sampler> {

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
      : MemUsageSamplerBase<SysInfo::Sampler>(name) {

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
