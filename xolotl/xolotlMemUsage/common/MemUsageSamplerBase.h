#ifndef XMEMUSAGE_MEM_USAGE_SAMPLER_BASE_H
#define XMEMUSAGE_MEM_USAGE_SAMPLER_BASE_H

#include "xolotlCore/Identifiable.h"
#include "xolotlMemUsage/IMemUsageSampler.h"


namespace xolotlMemUsage {

template<class T>
class MemUsageSamplerBase : public IMemUsageSampler,
                            public T,
                            public xolotlCore::Identifiable {

public:

    /// Disallow construction of a MemUsageSamplerBase without a name.
    MemUsageSamplerBase(void) = delete;

	/**
     * Construct a MemUsageSampler with the given name.
	 *
	 * @param name The object's name.
	 */
	MemUsageSamplerBase(const std::string& name)
      : xolotlCore::Identifiable(name),
        T(name) {

        // Nothing else to do.
    }

	/**
	 * Destroy the MemUsageSampler.
	 */
	virtual ~MemUsageSamplerBase(void) {

        // Ensure we've stopped sampling.
        stop();
    }

	/**
	 * Start sampling.
	 */
	virtual void start(void) {

        T::StartSampling();
    }

	/**
	 * Stop sampling.
	 */
	virtual void stop(void) {

        T::StopSampling();
    }
};

} // end namespace xolotlMemUsage

#endif // XMEMUSAGE_MEM_USAGE_SAMPLER_BASE_H
