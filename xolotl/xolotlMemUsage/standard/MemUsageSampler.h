#ifndef XMEMUSAGE_MEM_USAGE_SAMPLER_H
#define XMEMUSAGE_MEM_USAGE_SAMPLER_H

#include "xolotlMemUsage/common/MemUsageSamplerBase.h"
#include "xolotlMemUsage/standard/StatmSampler.h"


namespace xolotlMemUsage {

class MemUsageSampler : public MemUsageSamplerBase<Statm::Sampler> {

public:
    /**
     * Disallow building an object without a name.
     */
    MemUsageSampler(void) = delete;


	/**
     * Construct a MemUsageSampler with the given name.
	 *
	 * @param name The object's name.
	 */
	MemUsageSampler(const std::string& name)
      : MemUsageSamplerBase<Statm::Sampler>(name) {

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

#endif // XMEMUSAGE_MEM_USAGE_SAMPLER_H
