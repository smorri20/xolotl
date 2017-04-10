#ifndef XMEMUSAGE_MEM_USAGE_SAMPLER_H
#define XMEMUSAGE_MEM_USAGE_SAMPLER_H

#include "xolotlCore/Identifiable.h"
#include "xolotlMemUsage/IMemUsageSampler.h"
#include "xolotlMemUsage/standard/StatmSampler.h"


namespace xolotlMemUsage {

class MemUsageSampler : public IMemUsageSampler,
                            public StatmSampler,
                            public xolotlCore::Identifiable {

private:
    /**
     * Construct a MemUsageSampler.
     * Declared private to enforce that MemUsageSamplers must
     * be constructed with a name.
     */
    MemUsageSampler(void)
      : xolotlCore::Identifiable("unused"),
        StatmSampler("unused") {
        // Nothing else to do.
    }
    
public:

	/**
     * Construct a MemUsageSampler with the given name.
	 *
	 * @param name The object's name.
	 */
	MemUsageSampler(const std::string& name)
      : xolotlCore::Identifiable(name),
        StatmSampler(name) {

        // Nothing else to do.
    }

	/**
	 * Destroy the MemUsageSampler.
	 */
	virtual ~MemUsageSampler(void) {

        // Ensure we've stopped sampling.
        stop();
    }

	/**
	 * Start sampling.
	 */
	virtual void start(void) {

        StartSampling();
    }

	/**
	 * Stop sampling.
	 */
	virtual void stop(void) {

        StopSampling();
    }

	/**
	 * Obtain the sampler's metric values.
	 */
	virtual IMemUsageSampler::ValType getValue(void) const {

        return GetRunningSampleData().GetCurrentStats();
    }
};

} // end namespace xolotlMemUsage

#endif // XMEMUSAGE_MEM_USAGE_SAMPLER_H
