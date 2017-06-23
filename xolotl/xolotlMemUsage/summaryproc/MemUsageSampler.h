#ifndef XMEMUSAGE_MEM_USAGE_SAMPLER_H
#define XMEMUSAGE_MEM_USAGE_SAMPLER_H

#include "xolotlMemUsage/memUsageConfig.h"
#include "xolotlMemUsage/common/MemUsageSamplerBase.h"

#if defined(HAVE_STATM)
    #include "xolotlMemUsage/summaryproc/Statm/StatmSampler.h"
    namespace PerProcDataSource = xolotlMemUsage::Statm;
#else
    #error "Configuration error: thought we had a per-proc data source, but no recognized source available."
#endif // defined(HAVE_STATM)


namespace xolotlMemUsage {

class MemUsageSampler : public MemUsageSamplerBase<PerProcDataSource::Sampler> {

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
      : MemUsageSamplerBase<PerProcDataSource::Sampler>(name) {

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
